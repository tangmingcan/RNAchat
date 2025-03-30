from functools import lru_cache

import random
import string
from .constants import GeneID_URL, NUMBER_CPU_LIMITS
from .models import (
    CustomUser,
    userData,
    MetaFileColumn,
    SharedFile,
    UploadedFile,
)
from harmony import harmonize
import pandas as pd
import numpy as np
import scanpy as sc
import collections
from combat.pycombat import pycombat
from matplotlib import pyplot as plt
from django.http import HttpResponse, JsonResponse
from collections import Counter
import requests
import json
import io
import base64
from threadpoolctl import threadpool_limits
from sklearn.preprocessing import MinMaxScaler
from django.shortcuts import render
from functools import wraps
import jwt
import datetime
from django.conf import settings
from django.shortcuts import redirect
from datetime import timedelta
from django.utils.timezone import now
from django.contrib.auth import get_user_model
from sica.base import StabilizedICA
from rnanorm import CPM
from scipy.stats import gmean
import re
import ast


def delete_inactive_accounts():
    User = get_user_model()
    threshold_date = now() - timedelta(days=90)
    inactive_users = User.objects.filter(last_login__lt=threshold_date)
    count = inactive_users.count()
    inactive_users.delete()


# from pydeseq2.dds import DeseqDataSet
# from pydeseq2.default_inference import DefaultInference
# inference = DefaultInference(n_cpus=NUMBER_CPU_LIMITS)


def auth_required(f):
    @wraps(f)
    def wrap(request, *args, **kwargs):
        is_browser = "Mozilla" in request.META.get("HTTP_USER_AGENT", "")
        if is_browser:
            if not request.user.is_authenticated:
                return redirect("login")
        else:
            token = request.headers.get("Authorization")
            if not token or not token.startswith("Bearer "):
                return JsonResponse(
                    {"error": "Token is missing or invalid"}, status=401
                )
            token = token.split()[1]
            try:
                decoded_token = jwt.decode(
                    token, settings.SECRET_KEY, algorithms=["HS256"]
                )
                request.user = CustomUser.objects.get(
                    username=decoded_token["username"]
                )
            except (
                jwt.ExpiredSignatureError,
                jwt.InvalidTokenError,
                CustomUser.DoesNotExist,
            ):
                return JsonResponse({"error": "Invalid Token"}, status=401)
        return f(request, *args, **kwargs)

    return wrap


def generate_jwt_token(user):
    expiration = datetime.datetime.utcnow() + datetime.timedelta(hours=2)
    token = jwt.encode(
        {"username": user.username, "exp": expiration},
        settings.SECRET_KEY,
        algorithm="HS256",
    )
    return token


@lru_cache(maxsize=None)
def build_GOEnrichmentStudyNS():
    obo_fname = download_go_basic_obo()
    godag = GODag(obo_fname)
    gene2go_fname = download_ncbi_associations()
    gene2go_reader = Gene2GoReader(gene2go_fname, taxids=[9606])
    ns2assoc = gene2go_reader.get_ns2assc()  # bp,cc,mf
    return GOEnrichmentStudyNS(
        GENEID2NT.keys(),
        ns2assoc,
        godag,
        propagate_counts=False,
        alpha=0.05,
        methods=["fdr_bh"],
    )


def preview_dataframe(df):
    if len(df.columns) > 7:
        selected_columns = list(df.columns[:4]) + list(df.columns[-3:])
        df = df[selected_columns]
        df = df.rename(columns={df.columns[3]: "..."})
        df.loc[:, "..."] = "..."
    return df


def has_duplicate_columns(df):
    # Check for duplicated column names
    duplicated_columns = df.columns[df.columns.duplicated()]
    return len(duplicated_columns) > 0


def zip_for_vis(X2D, batch, obs):
    traces = {}
    for i, j, k in zip(X2D, batch, obs):
        j = str(j)
        if j not in traces:
            traces[j] = {"data": [i], "obs": [k]}
        else:
            traces[j]["data"].append(i)
            traces[j]["obs"].append(k)
    return traces


def fromPdtoSangkey(df):
    df = df.copy()
    source = []
    df[["parent_level"]] = [i + " " for i in df[["LABEL"]].values]
    columns = df.columns.tolist()

    for i in columns:
        source.extend([j[0] for j in df[[i]].values])
    nodes = collections.Counter(source)

    source = {}
    for i in range(1, len(columns)):
        agg = df.groupby([columns[i - 1], columns[i]]).size().reset_index()
        agg.columns = [columns[i - 1], columns[i], "count"]
        for index, row in agg.iterrows():
            source[(row[columns[i - 1]], row[columns[i]])] = row["count"]

    result = {}
    result1 = []
    dic_node = {}
    for i, j in enumerate(nodes.keys()):
        result1.append({"node": i, "name": j})
        dic_node[j] = i  # name->number
    result["nodes"] = result1

    result1 = []
    for i in source:
        result1.append(
            {"source": dic_node[i[0]], "target": dic_node[i[1]], "value": source[i]}
        )
    result["links"] = result1
    # return json.dumps(result)

    sou = []
    tar = []
    value = []
    result = {}
    for i in result1:
        sou.append(i["source"])
        tar.append(i["target"])
        value.append(i["value"])
    result["source1"] = sou
    result["target1"] = tar
    result["value1"] = value
    result["label"] = list(nodes.keys())
    return result


def go_it(test_genes):
    goeaobj = build_GOEnrichmentStudyNS()
    goea_results_all = goeaobj.run_study(
        [g.id for g in Gene.objects.filter(name__in=test_genes)]
    )
    goea_result_sig = [
        r for r in goea_results_all if r.p_fdr_bh < 0.05 and r.ratio_in_study[0] > 1
    ]
    go_df = pd.DataFrame(
        data=map(
            lambda x: [
                x.goterm.name,
                x.goterm.namespace,
                x.p_fdr_bh,
                x.ratio_in_study[0],
                GOTerm.objects.get(name=x.GO).gene.count(),
                list(GOTerm.objects.get(name=x.GO).gene.values_list("name", flat=True)),
            ],
            goea_result_sig,
        ),
        columns=["term", "class", "p_corr", "n_genes", "n_go", "genes"],
        index=map(lambda x: x.GO, goea_result_sig),
    )
    go_df["n_genes/n_go"] = go_df.n_genes / go_df.n_go
    return go_df


def get_random_string(length):
    letters = string.ascii_lowercase
    result_str = "".join(random.choice(letters) for i in range(length))
    return result_str


def combat(dfs):
    df_exp = pd.concat(dfs, join="inner", axis=1)
    batch = []
    datasets = dfs
    for j in range(len(datasets)):
        batch.extend([j for _ in range(len(datasets[j].columns))])
    df_corrected = pycombat(df_exp, batch)
    Xc = df_corrected.T
    return Xc.reset_index(drop=True)


# pip install harmony-pytorch
def harmony(dfs, batch, obs):
    dfs1 = pd.concat(dfs, join="inner", axis=0)
    adata = sc.AnnData(
        np.zeros(dfs1.values.shape), dtype=np.float64
    )  #'obs','FileName', 'LABEL'
    adata.X = dfs1.values
    adata.var_names = dfs1.columns.tolist()
    adata.obs_names = obs
    adata.obs["batch"] = batch
    # sc.tl.pca(adata)
    # sce.pp.harmony_integrate(adata, 'batch')
    # ha=adata.obsm['X_pca_harmony']
    # return pd.DataFrame(ha,columns=['PCA_'+str(i+1) for i in range(ha.shape[1])]).reset_index(drop=True)
    ha = harmonize(adata.X, adata.obs, batch_key="batch")
    return pd.DataFrame(ha, columns=dfs1.columns).reset_index(drop=True)


def bbknn(dfs):
    return dfs


def clusteringPostProcess(X2D, adata, method, usr):
    if method != "kmeans" and len(set(adata.obs[method])) == 1:
        # throw error for just 1 cluster
        raise Exception("Only 1 Cluster after clustering")

    li = adata.obs[method].tolist()
    count_dict = Counter(li)
    for member, count in count_dict.items():
        if count < 2:
            raise Exception(
                "The number of data in the cluster "
                + str(member)
                + " is less than 2, which will not be able for further analysis."
            )

    traces = zip_for_vis(X2D, list(adata.obs[method]), adata.obs_names.tolist())

    adata.obs["cluster"] = li
    usr.setAnndata(adata)

    barChart1 = []
    barChart2 = []

    with plt.rc_context():
        figure1 = io.BytesIO()
        figure2 = io.BytesIO()
        with threadpool_limits(limits=NUMBER_CPU_LIMITS, user_api="blas"):
            sc.tl.rank_genes_groups(adata, groupby=method, method="t-test")
            sc.tl.dendrogram(adata, groupby=method)
            sc.pl.rank_genes_groups_dotplot(
                adata, n_genes=4, show=False, color_map="bwr"
            )
            plt.savefig(figure1, format="svg", bbox_inches="tight")
            sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
            plt.savefig(figure2, format="svg", bbox_inches="tight")
    markers = sc.get.rank_genes_groups_df(adata, None)
    usr.setMarkers(markers)
    if usr.save() is False:
        raise Exception("Error for creating user records.")
    b = (
        adata.obs.sort_values(["batch1", method])
        .groupby(["batch1", method])
        .count()
        .reset_index()
    )

    b = b[["batch1", method, "batch2"]]
    b.columns = ["batch", method, "count"]
    barChart1 = [
        {
            "x": sorted(list(set(b[method].tolist()))),
            "y": b[b["batch"] == i]["count"].tolist(),
            "name": i,
            "type": "bar",
        }
        for i in set(b["batch"].tolist())
    ]

    b = (
        adata.obs.sort_values(["batch2", method])
        .groupby(["batch2", method])
        .count()
        .reset_index()
    )
    b = b[["batch2", method, "batch1"]]
    b.columns = ["batch", "cluster", "count"]
    barChart2 = [
        {
            "x": sorted(list(set(b["cluster"].tolist()))),
            "y": b[b["batch"] == i]["count"].tolist(),
            "name": i,
            "type": "bar",
        }
        for i in set(b["batch"].tolist())
    ]
    return {
        "traces": traces,
        "fileName": base64.b64encode(figure1.getvalue()).decode("utf-8"),
        "fileName1": base64.b64encode(figure2.getvalue()).decode("utf-8"),
        "bc1": barChart1,
        "bc2": barChart2,
        "method": usr.redMethod,
    }


def getTopGeneCSV(adata, groupby, n_genes):
    if len(set(adata.obs[groupby])) > 1:
        with threadpool_limits(limits=NUMBER_CPU_LIMITS, user_api="blas"):
            sc.tl.rank_genes_groups(adata, groupby=groupby, method="t-test")

        result = adata.uns["rank_genes_groups"]
        groups = result["names"].dtype.names  # Get the cluster names
        all_info = []

        # Extract all relevant information for each group
        for group in groups:
            filtered_genes = pd.DataFrame(
                {
                    f"Gene_{group}": result["names"][group],
                    f"LogFoldChange_{group}": result["logfoldchanges"][group],
                    f"PValue_{group}": result["pvals"][group],
                    f"AdjPValue_{group}": result["pvals_adj"][group],
                }
            )
            filtered_genes[f"Cluster_{group}"] = (
                group  # Add the cluster name for reference
            )

            # Sort by AdjPValue first (ascending order)
            filtered_genes = filtered_genes.sort_values(
                by=f"AdjPValue_{group}", ascending=True
            ).reset_index(drop=True)

            # Add the sorted DataFrame to the list
            all_info.append(filtered_genes)

        # Combine all clusters' data into one DataFrame horizontally
        result_df = pd.concat(all_info, axis=1)

        # Reset the index to ensure sequential indexing
        result_df = result_df.reset_index(drop=True)

        return result_df
    else:
        raise Exception("Only one cluster.")


def GeneID2SymID(geneList):
    geneIDs = []
    for g in geneList:
        if g.startswith("ENSG"):
            geneIDs.append(g)
    if len(geneIDs) == 0:
        return None
    ids_json = json.dumps(geneIDs)
    body = {"api": 1, "ids": ids_json}
    try:
        response = requests.post(GeneID_URL, data=body)
        text = json.loads(response.text)
    except:
        return None
    retRes = []
    for i in geneList:
        if i in text and text[i] is not None:
            retRes.append(text[i])
        else:
            retRes.append(i)
    return retRes


def UploadFileColumnCheck(df, keywords):
    for keyword in keywords:
        if keyword not in df.columns:
            return False, f"No {keyword} column in the expression file"
    if has_duplicate_columns(df):
        return False, "file has duplicated columns."
    return True, ""


def usrCheck(request, flag=1):
    username = request.user.username
    clientID = request.GET.get("cID", None)
    if clientID is None:
        return {"status": 0, "message": "clientID is Required."}
    if flag == 0:
        usr = userData(clientID, username)
    else:
        usr = userData.read(username, clientID)
        if usr is None:
            return {
                "status": 0,
                "message": "Can't find the user/device.Please request from the beginning.",
            }
    return {"status": 1, "usrData": usr}




# normalize clinic data
def normalize1(df, log2="No"):
    # Identify numeric columns
    return df

def loadSharedData(request, cID):
    files = UploadedFile.objects.filter(user=request.user, type1="exp", cID=cID).all()
    if files:
        files = [i.file.path for i in files]
    return files

def integrateCliData(request, cID):
    f = UploadedFile.objects.filter(user=request.user, type1="cli", cID=cID).first()
    if f is not None:
        temp0 = pd.read_csv(f.file.path)
        temp0 = temp0.dropna(axis=1)
    else:
        temp0 = pd.DataFrame()
    return temp0


def integrateExData(files, temp0):
    dfs = []
    batch = []
    obs = []

    for file in files:
        temp01 = pd.read_csv(file).set_index("ID_REF")
        # filter records that have corresponding clinical data in meta files.
        temp = temp01.join(temp0[["ID_REF", "LABEL"]].set_index("ID_REF"), how="inner")
        temp1 = temp.reset_index().drop(
            ["ID_REF", "LABEL"], axis=1, inplace=False
        )  # exclude ID_REF & LABEL
        temp1.index = temp.index
        # exclude NA
        temp1.dropna(axis=1, inplace=True)
        if temp1.shape[0] != 0:
            dfs.append(temp1)
            # color2.extend(list(temp.LABEL))
            batch.extend(
                ["_".join(file.split("_")[1:]).split(".csv")[0]]
                * temp1.shape[0]
            )
            # obs.extend(temp.ID_REF.tolist())
            obs.extend(temp1.index.tolist())
    

    dfs1 = dfs[0]
    dfs1["ID_REF"] = obs
    dfs1["FileName"] = batch
    return dfs1


def center_df(df):
    return df.apply(lambda x: x - x.mean(), axis=0)

"""flag==1 means response after uploading; flag==0 means initialize for the tab view."""


def getExpression(request, cID, flag=1):
    context = {}
    for file in UploadedFile.objects.filter(
        user=request.user, type1="exp", cID=cID
    ).all():
        df = pd.read_csv(file.file.path, nrows=5, header=0)
        f = file.file.name.split("/")[-1]

        check, mess = UploadFileColumnCheck(df, ("ID_REF",))
        if check is False:
            return HttpResponse(mess, status=400)
        df = preview_dataframe(df)
        json_re = df.reset_index().to_json(orient="records")
        data = json.loads(json_re)
        context["_".join(f.split("_")[1:])] = {}
        context["_".join(f.split("_")[1:])]["d"] = data
        context["_".join(f.split("_")[1:])]["names"] = ["index"] + df.columns.to_list()
    if flag == 0 and len(context) == 0:
        return HttpResponse("", status=200)
    return render(request, "table.html", {"root": context})


def getMeta(request, cID, flag=1):
    context = {}
    f = UploadedFile.objects.filter(user=request.user, type1="cli", cID=cID).first()
    if f is None:
        if flag == 1:
            return HttpResponse(mess, status=400)
        else:
            return HttpResponse("", status=200)
    context["metaFile"] = {}
    df = pd.read_csv(f.file.path, nrows=5, header=0)
    check, mess = UploadFileColumnCheck(df, ("ID_REF", "LABEL"))
    if check is False:
        return HttpResponse(mess, status=400)
    df = preview_dataframe(df)
    json_re = df.reset_index().to_json(orient="records")
    data = json.loads(json_re)
    context["metaFile"]["d"] = data
    context["metaFile"]["names"] = ["index"] + df.columns.to_list()
    return render(request, "table.html", {"root": context})


from sklearn.preprocessing import StandardScaler


def ICA(df11, df12=None, num=7, suffix=''):
    ta11 = df11.copy()
    sICA = StabilizedICA(n_components=num, n_runs=30, plot=False, n_jobs=-1)
    # ta11 = ta11.apply(func=lambda x: x - x.mean(), axis=0)
    sICA.fit(ta11)
    metageneCompose = pd.DataFrame(
        sICA.S_,
        columns=ta11.columns.values,
        index=["metapathway_" + str(i) for i in range(sICA.S_.shape[0])],
    )
    metagenes = sICA.transform(ta11)
    metagenes = pd.DataFrame(
        metagenes, columns=[i for i in metageneCompose.index], index=ta11.index
    )
    if df12 is None:
        return metageneCompose, metagenes, None, sICA
    metagenes2 = sICA.transform(df12)
    metagenes2 = pd.DataFrame(
        metagenes2, columns=[i for i in metageneCompose.index], index=df12.index
    )
    return metageneCompose, metagenes, metagenes2, sICA
    

def getColumnFields(numeric, request, usr):
    if numeric is None:
        result = [
            [i[0], i[1], i[2]]
            for i in MetaFileColumn.objects.filter(user=request.user, cID=usr.cID)
            .all()
            .values_list("colName", "label", "numeric")
            ]
    else:
        result = [
            [i[0], i[1], i[2]]
            for i in MetaFileColumn.objects.filter(
                user=request.user, cID=usr.cID, numeric=numeric
            )
            .all()
            .values_list("colName", "label", "numeric")
            ]
    return result
    

def expression_clinic_split(adata):
    df = adata.to_df()
    dfe = df.loc[:, ~df.columns.str.startswith("c_")]
    dfc = df.loc[:, df.columns.str.startswith("c_")]
    return dfe, dfc
    

def normalize_string(s):
    # Use regex to remove non-alphabetic characters and convert to lowercase
    return re.sub(r'[^a-zA-Z_]', '', s).lower()

def map_normalized_keys(keys, values):
    # Normalize all values and map them to their original
    normalized_map = {normalize_string(v): v for v in values}
    # Map each key to its corresponding original value if found
    return {k: normalized_map.get(normalize_string(k), None) for k in keys}

def drop_parallel(temp3, threshold=0.65):
    corr = temp3.corr()
    high_corr_pairs = []
    for i in range(len(corr.columns)//2):
        for j in range(i+len(corr.columns)//2,len(corr.columns)):
            if abs(corr.iloc[i, j]) > threshold:
                high_corr_pairs.append((corr.columns[i], corr.columns[j]))
    to_drop = set()
    for pair in high_corr_pairs:
        to_drop.add(pair[1])
    temp3_reduced = temp3.drop(columns=to_drop)
    return temp3_reduced
    
    
def MLparamSetting(mlMethod, ml_params):
    # Remove whitespace and split parameters
    cleaned_params = ml_params.replace(" ", "")
    param_pairs = cleaned_params.split(',')
    params_dict = {}
    
    # Parse parameters safely using ast.literal_eval
    if ml_params!='':
        try:
            for pair in param_pairs:
                key, value = pair.split('=')
                if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
                    params_dict[key] = value[1:-1]  # Remove quotes around the string
                else:
                    params_dict[key] = ast.literal_eval(value)  # Safely evaluate if it's not a string
        except (ValueError, SyntaxError):
            raise Exception(f"Invalid parameter format: '{pair}'")

    # Define default values
    defaults = {
        "RF": {"random_state":42},
        "XGB": {"random_state":42},
        "LIGHTGBM": {"random_state":42}
    }

    # Determine model-specific defaults
    if mlMethod == 'rf':
        model_defaults = defaults["RF"]
    elif mlMethod == 'xgb':
        model_defaults = defaults["XGB"]
    else:
        model_defaults = defaults["LIGHTGBM"]

    # Combine defaults with user-provided parameters
    params = {**model_defaults, **params_dict}

    return params

def shared_ml_dml_cf(usr, logY, mlType, label, cla, Tlabel='', Tclass=''):
    adata = usr.getAnndata()  
    temp = adata.to_df()
    dfe, dfc = expression_clinic_split(adata)
    if label=='LABEL':
        label='batch2'
    if Tlabel=='LABEL':
        Tlabel='batch2'
    
    if len(usr.metagenes)==0:
        raise Exception("Run ICA first.")
    X = usr.metagenes[0]
    if Tlabel=='':
        T = None
    else:
        T = adata.obs.loc[X.index, Tlabel].to_list()
        T = [1 if value == Tclass else 0 for value in T]

    if logY=='yes' and mlType!='classification':
        dfc = dfc.drop(columns=[label])
    if 'cell_type' in dfc.columns:
        dfc = dfc.drop(columns=['cell_type'])
    X = pd.concat([X, dfc],axis=1, join='inner')
 
    if mlType=='classification':
        if cla is None or label not in adata.obs_keys():
            raise Exception("incorrect message request")
        if cla =='cell_type':
            raise Exception("cell types can't be used for classification")
        dics = map_normalized_keys([cla], set(adata.obs[label].to_list()))
        target_value = dics.get(cla)
        if target_value is None:
            raise Exception(f"No such value: '{cla}' for column: '{label}'")
        y = adata.obs.loc[X.index, label].to_list()
        y = [1 if value == target_value else 0 for value in y]

    else:
        if logY is None or label not in temp.columns:
            raise Exception("incorrect message request")    
        y = temp.loc[X.index, label].to_list()
        if logY=='yes':
            y = np.log1p(y)
    return X, y, T
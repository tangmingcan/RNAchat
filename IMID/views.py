from django.shortcuts import render
from django.http import HttpResponse, JsonResponse

# from django.contrib.auth.decorators import login_required
import pandas as pd
import json
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import umap.umap_ as umap
import scanpy as sc
import numpy as np
from matplotlib import pyplot as plt, cm
import hdbscan
from threadpoolctl import threadpool_limits
from sklearn.preprocessing import StandardScaler
import math
import io
import base64
from django.db.models import Q
from collections import defaultdict
from django.contrib.auth import authenticate
import seaborn as sns
from pyecharts import options as opts
from pyecharts.charts import Graph

# gene_result.txt, genes_ncbi_proteincoding.py, go-basic.obo

import matplotlib
import re, json, shap
from sklearn.metrics import roc_curve, roc_auc_score, mean_squared_error

matplotlib.use("agg")
import plotly.graph_objects as go

import requests
from bs4 import BeautifulSoup
from .constants import BUILT_IN_LABELS, NUMBER_CPU_LIMITS, ALLOW_UPLOAD, SHAP_PLOT_DATASET
from .utils import (
    zip_for_vis,
    fromPdtoSangkey,
    GeneID2SymID,
    usrCheck,
    getExpression,
    getMeta,
    generate_jwt_token,
    auth_required,
    getColumnFields,
    expression_clinic_split,
    map_normalized_keys,
    MLparamSetting,
    shared_ml_dml_cf,
)

from .models import MetaFileColumn, UploadedFile, SharedFile, Msigdb
from django.db import transaction
from sklearn.model_selection import train_test_split
from IMID.tasks import (
    runIntegrate,
    runDgea,
    runClustering,
    runFeRed,
    runICA,
    runPCA,
    runTopFun,
    calculate_imf_scores,
    run_ora_task,
    run_oraPlot,
    run_Rf,
    run_XGBoost,
    run_LightGBM,
    run_shap,
    run_MLpreprocess,
    run_DMLreport,
    run_DMLreportDis,
    run_CFreport1,
    run_CFreport2,
    run_GB,
    run_matrixplot
)

from econml.dml import CausalForestDML

def restLogin(request):
    if request.method != "POST":
        return HttpResponse("Method not allowed.", status=405)
    body = request.body.decode("utf-8")
    try:
        data = json.loads(body)
    except:
        return JsonResponse({"error": "Invalid credentials"}, status=400)
    username = data["username"]
    password = data["password"]
    user = authenticate(username=username, password=password)
    if user is not None:
        token = generate_jwt_token(user)
        return JsonResponse({"token": token})
    return JsonResponse({"error": "Invalid credentials"}, status=400)


@auth_required
def index(request):
    return render(
        request,
        "index.html",
    )


"""
This page will show different shared data to different users within different user groups. If the user is defined in a group,
the programme will filter his corresponding group. Otherwise it will return all shared data.
"""


@auth_required
def tab(request):
    context = {}
    if request.user.groups.exists():
        context["cohorts"] = list(
            SharedFile.objects.filter(
                type1="expression", groups__in=request.user.groups.all()
            )
            .all()
            .values_list("cohort", "label")
        )
    else:
        context["cohorts"] = list(
            SharedFile.objects.filter(type1="expression")
            .all()
            .values_list("cohort", "label")
        )
    return render(request, "tab1.html", {"root": context})
      

@auth_required        
def compICAcohort(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    if usr.ica_cohort==():
        return HttpResponse('Please Run ICA based cohort first!', status=400)
    return JsonResponse({"text": usr.ica_cohort[:3]})   
        
@auth_required        
def compICAall(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    if usr.ica_cohort==():
        ica_cohort=[]
        X11 = usr.metagenes[0].columns.to_list()
        X22 = ''
    else:
        X11, X22 =usr.metagenes[0].columns.to_list(), usr.metagenes[3].columns.to_list()
        ica_cohort = usr.ica_cohort[:3]
    #temp1=pd.concat([X11,X22],axis=1)#0
    return JsonResponse({"ica_cohort": ica_cohort,"ica_cohort_metagenes":(X11, X22)})  



"""
This program provides entrance for user to upload/get Expression files. One user can upload multiple expression files at a time.
'ID_REF' should be one of the colnames in each Expression file.
User uses UploadedFile to store their uploaded expression data, the column type1 is set to be 'exp' to represent the stored expression file.
"""


@auth_required
def opExpression(request):
    if request.method == "POST":
        if ALLOW_UPLOAD is False:
            return HttpResponse("Not Allowed user upload", status=400)
        cID = request.POST.get("cID", None)
        if cID is None:
            return HttpResponse("cID not provided.", status=400)
        UploadedFile.objects.filter(user=request.user, type1="exp", cID=cID).delete()
        files = request.FILES.getlist("files[]", None)
        new_exp_files = []

        f=files[0]
        temp_file = UploadedFile(user=request.user, cID=cID, type1="exp", file=f)
        new_exp_files.append(temp_file)
        UploadedFile.objects.bulk_create(new_exp_files)
        return getExpression(request, cID)
    elif request.method == "GET":
        cID = request.GET.get("cID", None)
        if cID is None:
            return HttpResponse("cID not provided.", status=400)
        return getExpression(request, cID, 0)


"""
This program provides entrance for user to upload/get Meta files. One user can upload only one meta file at a time.
'ID_REF','LABEL' should be one of the colnames in each Meta file.
User uses UploadedFile to store their uploaded expression data, the column type1 is set to be 'cli' to represent the stored clinical data.
"""


@auth_required
def opMeta(request):
    if request.method == "POST":
        if ALLOW_UPLOAD is False:
            return HttpResponse("Not Allowed user upload", status=400)
        cID = request.POST.get("cID", None)
        if cID is None:
            return HttpResponse("cID not provided.", status=400)
        files = request.FILES.getlist("meta", None)
        if files is None:
            return HttpResponse("Upload the meta file is required", status=405)
        files = files[0]
        UploadedFile.objects.filter(user=request.user, type1="cli", cID=cID).delete()
        UploadedFile.objects.create(user=request.user, cID=cID, type1="cli", file=files)
        return getMeta(request, cID)
    elif request.method == "GET":
        cID = request.GET.get("cID", None)
        if cID is None:
            return HttpResponse("cID not provided.", status=400)
        return getMeta(request, cID, 0)


"""
This is for data integration based on selected options. The files comes from 1. user uploaded files. 2. Built-in Data.
in_ta and in_ta1 are used for expression and meta data list separately
The processing logic is:
First, using user uploaded meta file as the base then inner join with selected built-in meta files to get the shared clinic features. =>temp0
Then join the expression files=>dfs1 with consideration of log2, batch effect
Third, join temp0 and dfs1;
"""


@auth_required
def downloadICA(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    typ = request.GET.get("type", 'ica')
    if len(usr.metagenes)!=1:
        return HttpResponse("Please run ICA first.", status=400)
    result_df = usr.metagenes[0]
    response = HttpResponse(content_type="text/csv")
    response["Content-Disposition"] = "attachment; filename=ica.csv"
    result_df.to_csv(path_or_buf=response)
    return response



@auth_required
def checkUser(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        return HttpResponse("User exists.", status=200)


@auth_required
def meta_columns(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    cID = request.GET.get("cID", None)
    if request.method == "GET":
        numeric = request.GET.get("numeric", None)
        param = request.GET.get("param",1)
        if numeric is None:
            result = [
                [i[0], i[1], i[2]]
                for i in MetaFileColumn.objects.filter(user=request.user, cID=usr.cID)
                .all()
                .values_list("colName", "label", "numeric") if i[0]!='cell_type'
            ]
            if param=="2":
                created =[i[0].split('__crted')[0] for i in result if 'crted' in i[0]]
                result = [i for i in result if i[0] not in created]
        else:
            result = [
                [i[0], i[1], i[2]]
                for i in MetaFileColumn.objects.filter(user=request.user, cID=usr.cID, numeric="0")
                .all().values_list("colName", "label", "numeric") if i[0]!='cell_type'
            ]
            if param =="2":
                created =[i[0].split('__crted')[0] for i in result if 'crted' in i[0]]
                result2 = [
                    [i[0], i[1], i[2]]
                    for i in MetaFileColumn.objects.filter(
                        user=request.user, cID=usr.cID, numeric=numeric
                    )
                    .all()
                    .values_list("colName", "label", "numeric")
                ]       
                result = [i for i in result2 if i[0] not in created]
        return JsonResponse(result, safe=False)
    elif request.method == "PUT":
        labels = request.GET.get("labels", None)
        if labels is None:
            return HttpResponse("Labels illegal.", status=400)
        labels = [i.strip() for i in labels.split(",")]
        error_labels = BUILT_IN_LABELS.intersection(set(labels))
        if len(error_labels) != 0:
            return HttpResponse(
                "Labels creating Problem, can't use retained words as an label:"
                + str(error_labels),
                status=400,
            )
        try:
            with transaction.atomic():
                df = usr.getIntegrationData().copy()
                crtedDic = defaultdict(list)
                for i in df.columns:
                    if "__crted" in i:
                        crtedDic[i.split("__crted")[0] + "__crted"].append(i)
                if labels != [""]:
                    labels_t = labels.copy()
                    labels1 = set(
                        [
                            i.split("__crted")[0] + "__crted"
                            for i in labels
                            if "__crted" in i
                        ]
                    )
                    labels2 = [i for i in labels if "__crted" not in i]
                    labels = set()
                    for i in crtedDic:
                        if i in labels1:
                            labels.update(crtedDic[i])  # add age__crted1-N
                            labels.add(i.split("__crted")[0])  # add age
                        else:
                            labels.update(crtedDic[i])  # add agg_crted1-N
                    labels.update(labels2)
                    df1 = df.drop(labels, axis=1, inplace=False)
                    usr.setAnndata(df1)
                    adata = usr.getAnndata()
                    for label in labels_t:
                        if not np.issubdtype(df[label].dtype, np.number):
                            adata.obs[label] = df[label]
                        else:
                            raise Exception(
                                "Can't auto convert numerical value for label."
                            )

                    MetaFileColumn.objects.exclude(colName="LABEL").exclude(
                        user=request.user, cID=cID, colName__in=labels_t
                    ).update(label="0")
                    MetaFileColumn.objects.filter(
                        user=request.user, cID=cID, colName__in=labels_t
                    ).update(label="1")
                else:
                    usr.setIntegrationData(df)
                    MetaFileColumn.objects.exclude(colName="LABEL").filter(
                        user=request.user, cID=cID
                    ).update(label="0")
        except Exception as e:
            return HttpResponse("Labels creating Problem. " + str(e), status=400)
        finally:
            if usr.save() is False:
                return HttpResponse("Can't save user record", status=500)
        return HttpResponse("Labels updated successfully.", status=200)
    elif request.method == "POST":
        post_data = json.loads(request.body)
        colName = post_data.get("colName")
        threshold = post_data.get("thredshold")
        threshold_labels = post_data.get("thredshold_labels")
        if (
            colName is None
            or threshold is None
            or threshold_labels is None
            or len(threshold) != len(threshold_labels) - 1
        ):
            return HttpResponse("Illegal Param Input. ", status=400)

        threshold = [float(i) for i in threshold]
        count = MetaFileColumn.objects.filter(
            (Q(colName=colName) | Q(colName__startswith=colName + "__crted"))
            & Q(user=request.user)
        ).count()
        if count == 0:
            HttpResponse(
                "No such colName: " + str(colName) + " in meta file.", status=400
            )
        colName1 = colName + "__crted" + str(count)
        df, adata = usr.getIntegrationData(), usr.getAnndata()
        conditions = [df[colName] < threshold[0]]
        for i in range(len(threshold_labels) - 2):
            conditions.append(
                (df[colName] >= threshold[i]) & (df[colName] < threshold[i + 1])
            )
        conditions.append(df[colName] >= threshold[-1])

        try:
            with transaction.atomic():
                df[colName1] = np.select(conditions, threshold_labels)
                adata.obs[colName1] = df[colName1].copy()
                MetaFileColumn.objects.create(
                    user=request.user, cID=cID, colName=colName1, label="0", numeric="0"
                )
        except Exception as e:
            return HttpResponse("Labels creating Problem. " + str(e), status=400)
        finally:
            if usr.save() is False:
                return HttpResponse("Can't save user record", status=500)
        return HttpResponse("Label created Successfully. ", status=200)

@auth_required
def meta_column_values(request, colName):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    adata = usr.getAnndata()
    df = usr.getIntegrationData()
    if colName.lower() == "Cluster".lower():
        colName = "cluster"
    if request.method == "GET":
        if colName == 'LABEL':
            colName = 'batch2'
        if colName in adata.obs_keys() and not np.issubdtype(
            adata.obs[colName].dtype, np.number
        ):
            temp = list(set(adata.obs[colName]))
            if len(temp) == 1:
                return HttpResponse(
                    "Only 1-type value found in the colName: " + colName, status=400
                )
            elif len(temp) > 30:
                return HttpResponse(
                    "More than 30-type values found in the colName: " + colName,
                    status=400,
                )
            temp.sort()
            return JsonResponse(temp, safe=False)
        if colName in df.columns:
            histogram_trace = go.Histogram(
                x=df[colName],
                histnorm="probability density",  # Set histogram normalization to density
                marker_color="rgba(0, 0, 255, 0.7)",  # Set marker color
            )

            # Configure the layout
            layout = go.Layout(
                title="Density Plot for " + colName,  # Set plot title
                xaxis=dict(title=colName),  # Set x-axis label
                yaxis=dict(title="Density"),  # Set y-axis label
            )

            # Create figure
            fig = go.Figure(data=[histogram_trace], layout=layout)

            return HttpResponse(
                base64.b64encode(fig.to_image(format="svg")).decode("utf-8"),
                content_type="image/svg+xml",
            )
        else:
            return HttpResponse("Can't find the colName: " + colName, status=400)
    if request.method == "DELETE":
        col = MetaFileColumn.objects.filter(
            user=request.user, cID=usr.cID, colName=colName
        ).first()
        if col is None:
            return HttpResponse("No such colName called:" + colName, status=400)
        elif col.label == "1":
            return HttpResponse("Please make {colName} inactive first.", status=400)
            
        try:
            with transaction.atomic():
                MetaFileColumn.objects.filter(
                    user=request.user, cID=usr.cID, colName=colName, label="0"
                ).delete()
                usr.integrationData = usr.integrationData.drop(columns=[colName,])
        except Exception as e:
            return HttpResponse("Labels creating Problem. " + str(e), status=400)
        finally:
            if usr.save() is False:
                return HttpResponse("Can't save user record", status=500)
        return HttpResponse("Delete {colName} Successfully.", status=200)
        
@auth_required
def edaIntegrate(request):
    checkRes = usrCheck(request, 0)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    cID = request.GET.get("cID", None)

    try:
        result = runIntegrate.apply_async(
            (request, cID, usr), serializer="pickle"
        ).get()
    except Exception as e:
        return HttpResponse(str(e), status=500)
    return HttpResponse("Operation successful.", status=200)

from urllib.parse import unquote
@auth_required
def ICAreport(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    cID = request.GET.get("cID", None)

    method = request.GET.get("method", 'ica')
    adata = usr.getAnndata().copy()
    dfe, dfc = expression_clinic_split(adata)
    try:
        dim = int(request.GET.get("dim",10))
    except Exception:
        return HttpResponse("Dim should be an integer", status=500)
    if dim >= dfe.shape[0]:
        return HttpResponse("Dimension is too large!" , status=500)     
    #if method=='ica':
    try:
        metageneCompose, metagenes, _, _ = runICA.apply_async((dfe, dim,''),serializer="pickle").get()
    except Exception as e:
        return HttpResponse("ICA Failed:" + str(e), status=500)
    usr.metagenes = [metagenes,]
    usr.metageneCompose = [metageneCompose,]
    usr.imf = None
    usr.ica_cohort = ()
    #print(usr.getAnndata().to_df().columns)
    #print(usr.getAnndata().obs.columns)
    #print(usr.getIntegrationData().columns)
    if usr.save() is False:
        return HttpResponse("Can't save user record", status=500)
    return HttpResponse("ICA Successful.", status=200)


@auth_required
def goML(request, mlType, mlMethod):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    #ml/classification/rf/?cID=12345&label=LABEL&class=nonrespondershareORBIT&testRatio=2
    #ml/regression/rf/?cID=12345&label=cage&logY=no&testRatio=2
    # print(len(usr.metagenes))
    try:
        body = json.loads(request.body)
        ml_param = body.get('ml_params', '')
        label = request.GET.get("label", None)
        cla = request.GET.get("class", None)
        testRatio=request.GET.get('testRatio',None)
        logY = request.GET.get('logY',None)
        mlType, param, X, y, testRatio, _ = run_MLpreprocess.apply_async((label, cla, testRatio, logY, usr, mlMethod, mlType,  ml_param),serializer="pickle").get()
        # print(mlType)
        # print(param)
        # print(X)
        #print(X.shape)
        #print(X)
        # print(y)
        # print(testRatio)
        # print(len(usr.metagenes))            
    except Exception as e:
        return HttpResponse(str(e), status=500)
    if mlType=='classification':    
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=testRatio, stratify=y, random_state=42)
        # print(X.shape,X_train.shape,X_test.shape)
    else:
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=testRatio, random_state=42)
    if mlMethod=='rf':
        result=run_Rf.apply_async((mlType, param, X, X_train, y_train, X_test, y_test,testRatio),serializer="pickle").get()
    elif mlMethod=='xgb':
        result=run_XGBoost.apply_async((mlType, param, X, X_train, y_train, X_test, y_test, testRatio),serializer="pickle").get()
    else:
        result=run_LightGBM.apply_async((mlType, param, X, X_train, y_train, X_test, y_test, testRatio),serializer="pickle").get()
    usr.model=[result[1],]
    if usr.save() is False:
        return HttpResponse("Can't save user record", status=500)
    return JsonResponse({"result": result[0]})
        
@auth_required
def goDML(request, mlType, mlMethod,number):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    #dml/classification/rf/1/?cID=12345&label=LABEL&class=nonrespondershareORBIT&testRatio=2
    #dml/regression/rf/1/?cID=12345&label=cage&logY=no&testRatio=2
    try:
        body = json.loads(request.body)
        ml_params = body.get('ml_params', '')
        label = request.GET.get("label", None)
        cla = request.GET.get("class", None)
        testRatio=request.GET.get('testRatio',None)
        logY = request.GET.get('logY',None)
        mlType, param, X, y, testRatio, _ = run_MLpreprocess.apply_async((label, cla, testRatio, logY, usr, mlMethod, mlType, ml_params),serializer="pickle").get()
        X['y'] =y
        temp = X.copy()
        adata = usr.getAnndata().copy()
        T = body['T']['batch']
        if T=='LABEL':
            T = 'batch2'
        if number==1:
            condition = adata.obs[T] == body['T']['batch1']
        else:
            condition = adata.obs[T] == body['T']['batch2']
        usr.ica_cohort = (T, body['T']['batch1'], body['T']['batch2'])
        valid_indices = condition[condition].index
        X = temp.loc[valid_indices].drop(columns=['y'])
        y = temp.loc[valid_indices,'y'].values
    except Exception as e:
        return HttpResponse(str(e), status=500)
    if mlType=='classification':    
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=testRatio, stratify=y, random_state=42)
    else:
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=testRatio, random_state=42)
    if mlMethod=='rf':
        result=run_Rf.apply_async((mlType, param, X, X_train, y_train, X_test, y_test, testRatio),serializer="pickle").get()
    elif mlMethod=='xgb':
        result=run_XGBoost.apply_async((mlType, param, X, X_train, y_train, X_test, y_test, testRatio),serializer="pickle").get()
    else:
        result=run_LightGBM.apply_async((mlType, param, X, X_train, y_train, X_test, y_test, testRatio),serializer="pickle").get()
    if number==1:
        usr.DMLmodel=[result[1],]
    else:
        usr.DMLmodel=[usr.DMLmodel[0],result[1]]
    if usr.save() is False:
        return HttpResponse("Can't save user record", status=500)
    return JsonResponse({"result": result[0]})

@auth_required
def feaImp(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    flip = request.GET.get("flip", 'no')
    if len(usr.model)==0:
        return HttpResponse("Please run ML first.", status=500)
    if len(usr.model)==1:
        image_data, usr.shapCache = run_shap.apply_async((usr.model[0][0], usr.model[0][1],flip),serializer="pickle").get()
    else:
        image_data = ''
    if usr.save() is False:
        return HttpResponse("Can't save user record", status=500)
    return render(request, "shap_importance.html", {"svg_content": image_data, "clientID":usr.cID})

@auth_required
def DMLfeaImp(request, number):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    if len(usr.DMLmodel)<number:
        return HttpResponse("Please run DML first.", status=500)
    if number==1:
        image_data, usr.shapCache = run_shap.apply_async((usr.DMLmodel[0][0], usr.DMLmodel[0][1]),serializer="pickle").get()
    else:
        image_data, usr.shapCache = run_shap.apply_async((usr.DMLmodel[1][0], usr.DMLmodel[1][1]),serializer="pickle").get()
    if usr.save() is False:
        return HttpResponse("Can't save user record", status=500)
    return render(request, "shap_importance.html", {"svg_content": image_data,"clientID":usr.cID})
        
@auth_required
def feaImpData(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    if usr.shapCache is None:
        return HttpResponse("Please run Shap Importance first.", status=500)
    method = request.GET.get("method", 'abs')
    shap_values = usr.shapCache
    if method =='abs':
        shap_importance = np.abs(shap_values.values).mean(axis=0)
    else:
        shap_importance = np.ptp(shap_values.values, axis=0)  # np.ptp computes peak-to-peak (max - min)
        
    feature_names = shap_values.feature_names
    # print(feature_names)
    result_df = pd.DataFrame({
        "Feature": feature_names,
        "SHAP_Importance": shap_importance
    })
    result_df.sort_values(by="SHAP_Importance", ascending=False, inplace=True)
    response = HttpResponse(content_type="text/csv")
    response["Content-Disposition"] = "attachment; filename=shap_importance.csv"
    result_df.to_csv(path_or_buf=response, index=False)
    return response
          
@auth_required
def DMLreport(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    if len(usr.DMLmodel)!=2:
        return HttpResponse("Please run 2 DML models first.", status=500)
    image_data, usr.shapImportance, shap_importance_top6 = run_DMLreport.apply_async((usr.DMLmodel[0][0], usr.DMLmodel[0][1], usr.DMLmodel[1][0],usr.DMLmodel[1][1]), serializer="pickle").get()
    if usr.save() is False:
        return HttpResponse("Can't save user record", status=500) 
    return render(request, "DMLshap_importance.html", {"svg_content": image_data, 
        "feature_names": ", ".join(shap_importance_top6),
        "clientID":usr.cID })
        
     
@auth_required
def DMLreportDis(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    param = request.GET.get("param", None)
    if param is None:
        return HttpResponse("No input Feature names", status=404)
    
    metagenes = param.replace("\t", "").replace(" ", "").split(',')
    # Assuming result_final1 contains one tuple with two SHAP importance sorted dictionaries
    if len(usr.shapImportance)!=2:
        return HttpResponse("Please run DML first.", status=500)
    shap_importance_sorted1 = pd.DataFrame(usr.shapImportance[0].values, columns=usr.DMLmodel[0][1].columns).abs().mean(axis=0)
    shap_importance_sorted2 = pd.DataFrame(usr.shapImportance[1].values, columns=usr.DMLmodel[1][1].columns).abs().mean(axis=0)
    
    image_data = run_DMLreportDis.apply_async((metagenes, shap_importance_sorted1, shap_importance_sorted2, usr.ica_cohort[:3]), serializer="pickle").get()    
    return JsonResponse({"image": image_data})
        
@auth_required
def ShapInteraction(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    else:
        usr = checkRes["usrData"]
    param = request.GET.get("param", None)
    metagenes = param.replace("\t", "").replace(" ", "").split(',')
    if len(metagenes)==0:
        return HttpResponse("No metagene", status=404)
    if len(metagenes)==1:
        if metagenes[0] not in usr.shapCache.feature_names:
            return HttpResponse("Illegal metagene", status=500)
        with threadpool_limits(limits=NUMBER_CPU_LIMITS, user_api="blas"):
            figure1 = io.BytesIO()
            shap.plots.scatter(usr.shapCache[:, metagenes[0]], color=usr.shapCache)
            plt.savefig(figure1, format="svg", bbox_inches="tight")
            plt.close()
    else:
        if metagenes[0] not in usr.shapCache.feature_names or metagenes[1] not in usr.shapCache.feature_names:
            return HttpResponse("Illegal metagene", status=500)
        with threadpool_limits(limits=NUMBER_CPU_LIMITS, user_api="blas"):
            figure1 = io.BytesIO()
            shap.plots.scatter(usr.shapCache[:, metagenes[0]], color=usr.shapCache[:, metagenes[1]])
            plt.savefig(figure1, format="svg", bbox_inches="tight")
            plt.close()
    image_data = base64.b64encode(figure1.getvalue()).decode("utf-8")
    return JsonResponse({"image": image_data})
        
@auth_required
def CFreport1(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)

    usr = checkRes["usrData"]

    try:
        params = json.loads(request.body)
        paramT = MLparamSetting(params.get('T_model', 'rf'), params.get('T_p', ''))
        paramY = MLparamSetting(params.get('Y_model', 'rf'), params.get('Y_p', ''))
        testRatio = float(params['testRatio'])
        X, y, T = shared_ml_dml_cf(usr, params['Y_log'], params['Y_type'], params['Y_label'], params['Y_interest'],  params['T_label'], params['T_interest'])
        #print(X)
        #print(y)
        #print(T)
    except Exception as e:
        return HttpResponse(f"Illegal parameters: {e}", status=500)
 
    mlType = params.get('Y_type', 'classification')

    X_train, X_test, T_train, T_test, y_train, y_test = train_test_split(X, T, y, test_size=testRatio, stratify=T, random_state=42)
    try:
        model_func_map = {'rf': run_Rf, 'xgb': run_XGBoost, 'lightgbm': run_LightGBM, 'gb':run_GB}
        run_treatment = model_func_map.get(params['T_model'], run_Rf)
        resultT = run_treatment.apply_async(('classification', paramT, X, X_train, T_train, X_test, T_test, testRatio), serializer="pickle").get()

        run_outcome = model_func_map.get(params['Y_model'], run_Rf)
        resultY = run_outcome.apply_async((mlType, paramY, X, X_train, y_train, X_test, y_test, testRatio), serializer="pickle").get()

        model_t, model_y = resultT[1][0], resultY[1][0]
        performance_image = run_CFreport1.apply_async((mlType, model_t, model_y, X_test, T_test, y_test), serializer="pickle").get()
    except Exception as e:
        return HttpResponse(f"Failed running CF: {e}", status=500)
        
    usr.model_ty = (model_t, model_y)
    if usr.save() is False:
        return HttpResponse("Can't save user record", status=500)
    return render(request, "imf_template.html", {"svg_content": performance_image, "title":"Model Performance", "message":""})

@auth_required
def CFreport2(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    usr = checkRes["usrData"]
    if len(usr.model_ty)==0:
        return HttpResponse('Run T&Y model First.', status=400)
    try:
        params = json.loads(request.body)
        testRatio = float(params['testRatio'])
        mlType = params.get('Y_type', 'classification')
        cf_params = params.get('cf_params', '')
        cf_params = MLparamSetting('rf', cf_params)
        model_t, model_y = usr.model_ty
        if mlType == 'classification':
            causal_forest = CausalForestDML(model_t=model_t, model_y=model_y, discrete_treatment=True, discrete_outcome=True, **cf_params)
        else:
            causal_forest = CausalForestDML(model_t=model_t, model_y=model_y, discrete_treatment=True, discrete_outcome=False, **cf_params)
            
        X, y, T = shared_ml_dml_cf(usr, params['Y_log'], params['Y_type'], params['Y_label'], params['Y_interest'], params['T_label'], params['T_interest'])
        # print(X)
        # print(y)
        # print(T)
        X_train, X_test, T_train, T_test, y_train, y_test = train_test_split(X, T, y, test_size=testRatio, stratify=T, random_state=42)
    except Exception as e:
        return HttpResponse(f"Error generating CFreport: {e}", status=500)

    causal_forest.fit(y_train, T_train, X=X_train)
    image = run_CFreport2.apply_async((causal_forest, X_test), serializer="pickle").get()
    usr.cf = (causal_forest, X, X_test.index.to_list())
    message = "Training ATE: "+str(causal_forest.ate_)
    message += "     ATE_STD_ERR: "+str(causal_forest.ate_stderr_)
    if usr.save() is False:
        return HttpResponse("Can't save user record", status=500)
    return render(request, "imf_template.html", {"svg_content": image, "title":"Treatment Result (Test dataset)", "message": message})
        
@auth_required  
def CFfeaImp(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    usr = checkRes["usrData"]
    if len(usr.cf)==0:
        return HttpResponse(f"Run Causal Forests first.", status=500)
    causal_forest = usr.cf[0]
    X, X_test = usr.cf[1], usr.cf[1].loc[usr.cf[2]]
    if SHAP_PLOT_DATASET=='X':
        shap_values = causal_forest.shap_values(X)
    elif SHAP_PLOT_DATASET=='X_train':
        shap_values = causal_forest.shap_values(X.drop(usr.cf[2]))
    else:
        shap_values = causal_forest.shap_values(X_test) 
    shap_val = shap_values['Y0']['T0_1']

    # Calculate spread (max - min) for each feature
    spread = np.ptp(shap_val.values, axis=0)  # np.ptp computes peak-to-peak (max - min)
    
    sorted_indices = np.argsort(spread)[::-1]
    figure = io.BytesIO()
    shap.plots.beeswarm(
        shap_val,
        max_display=35,
        show=False,
        order=sorted_indices
    )
    plt.savefig(figure, format="svg", bbox_inches="tight")
    plt.close()
    image = base64.b64encode(figure.getvalue()).decode("utf-8")
    usr.shapCache = shap_val
    if usr.save() is False:
        return HttpResponse("Can't save user record", status=500)
    return render(request, "shap_importance.html", {"svg_content": image,"clientID":usr.cID})

@auth_required
def dataVisual1(request):
    #usr.metagenes = [metagenes,]
    #usr.metageneCompose = [metageneCompose,]
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    usr = checkRes["usrData"]
    
    if len(usr.metageneCompose) == 0:
        return HttpResponse(f"Run ICA first.", status=500)
    
    thre = request.GET.get("thre", None)
    if thre is None:
        thre = 3
        dff = runTopFun.apply_async((usr.metageneCompose[0], thre), serializer="pickle").get()
        
        # Generate the correlation matrix
        corr_matrix = usr.metagenes[0].corr(method='pearson')
    
        # Plot the heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm', linewidths=0.5, cbar=True)
    
        # Save the plot to a buffer
        fig1 = io.BytesIO()
        plt.savefig(fig1, format='svg')
        fig1.seek(0)
        plot_svg = base64.b64encode(fig1.getvalue()).decode("utf-8")
        fig1.close()
        
        usr.ora = dff.reset_index()
        
        # Prepare the context to pass to the template
        context = {
            'root': {"df": dff.reset_index().to_dict(orient='records'), "cID": usr.cID},
            'plot_svg': plot_svg
        }
        if usr.save() is False:
            return HttpResponse("Can't save user record", status=500)
        return render(request, "dv1.html", context)
    else:
        try:
            thre = float(thre)
        except ValueError:
            return JsonResponse({"error": "Invalid 'thre' value, it must be a float."}, status=400)
        dff = runTopFun.apply_async((usr.metageneCompose[0], thre), serializer="pickle").get()   
        usr.ora = dff.reset_index()
        if usr.save() is False:
            return HttpResponse("Can't save user record", status=500)
        return JsonResponse(dff.reset_index().to_dict(orient='records'), safe=False)

@auth_required
def dataVisual2(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)
    
    usr = checkRes["usrData"]
    if usr.ora is None:
        return HttpResponse("Please run data Visualisation I first.", status=500)
    
    threshold = request.GET.get('threshold', None)
    cID = request.GET.get('cID', None)
    
    # If no threshold provided, render the HTML template
    if threshold is None:
        return render(request, 'dv2.html', {"cID": usr.cID})
    
    # Process with the threshold
    threshold = int(threshold)
    
    # Copy and prepare data
    df = usr.ora.copy()
    df.index.name = 'meta_id'
    df.reset_index(inplace=True)
    df_dict = df.to_dict(orient='records')
    df_dict_temp = {i['index']: i['inputs'] for i in df_dict}
    df_dict = df_dict_temp
    
    # Get unique metagenes and group them
    unique_metagenes = sorted([i for i in usr.metageneCompose[0].columns])
    group_order = sorted(set(['_' + i.split('_')[-1] for i in usr.metageneCompose[0].columns]))
    
    grouped_metagenes = {suffix: [] for suffix in group_order}
    for metagene in unique_metagenes:
        for suffix in group_order:
            if metagene.endswith(suffix):
                grouped_metagenes[suffix].append(metagene)
                break
    
    # Map metagenes to meta_ids
    metagene_to_meta = {}
    for meta_id, genes in df_dict.items():
        for gene in genes:
            if gene not in metagene_to_meta:
                metagene_to_meta[gene] = []
            metagene_to_meta[gene].append(meta_id)
    
    # Create edges and include shared meta groups information
    edges = []
    for i in range(len(unique_metagenes)):
        for j in range(i + 1, len(unique_metagenes)):
            meta_groups_1 = set(metagene_to_meta.get(unique_metagenes[i], []))
            meta_groups_2 = set(metagene_to_meta.get(unique_metagenes[j], []))
            shared_meta_groups = meta_groups_1.intersection(meta_groups_2)
            num_shared = len(shared_meta_groups)
            
            if num_shared >= threshold:
                edges.append({
                    'source': unique_metagenes[i],
                    'target': unique_metagenes[j],
                    'num_shared': num_shared,
                    'shared_meta_groups': shared_meta_groups  # Add shared_meta_groups information
                })
    # Create nodes and links with widths proportional to num_shared
    nodes = [{"name": m, "symbolSize": 10, "category": group_order.index('_' + m.split('_')[-1])} for m in unique_metagenes]
    links = [
        {
            "source": e['source'], 
            "target": e['target'], 
            "value": e['num_shared'],
            "lineStyle": {
                "width": e['num_shared'] * 2,  # Multiply by 2 to make width more visible
                "opacity": 0.7,
                "curveness": 0.3
            },
            "tooltip": {
                "formatter": f"Shared Meta Groups: {', '.join(e['shared_meta_groups'])} <br> Number of Shared: {e['num_shared']}"  # Display shared_meta_groups and num_shared
            }
        } for e in edges
    ]
    categories = [{"name": g} for g in group_order]
    
    # Create graph with pyecharts
    graph = (
        Graph(init_opts=opts.InitOpts(width="1000px", height="600px"))
        .add(
            "",
            nodes=nodes,
            links=links,
            categories=categories,
            layout="circular",
            is_rotate_label=True,
            linestyle_opts=opts.LineStyleOpts(curve=0.3),
            label_opts=opts.LabelOpts(position="right"),
        )
        .set_global_opts(
            title_opts=opts.TitleOpts(title="Metagene Network"),
            legend_opts=opts.LegendOpts(orient="vertical", pos_left="2%", pos_top="20%"),
            tooltip_opts=opts.TooltipOpts(
                formatter="{b}: {c}"  # The default tooltip formatter
            ),
        )
    )
    
    # Return the rendered graph HTML
    return JsonResponse({"graph_html": graph.render_embed()})

    
    
@auth_required
def matrixplot(request):
    checkRes = usrCheck(request)
    if checkRes["status"] == 0:
        return HttpResponse(checkRes["message"], status=400)

    usr = checkRes["usrData"]
    by = request.GET.get("by", 'batch2')
    temp = usr.metagenes[0].copy()
    image_base64 = run_matrixplot.apply_async((temp, usr, by), serializer="pickle").get()
    return JsonResponse({"image": image_base64})
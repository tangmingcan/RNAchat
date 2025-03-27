from .constants import NUMBER_CPU_LIMITS, ONTOLOGY, SHAP_PLOT_DATASET
import scanpy as sc
import math
from matplotlib import pyplot as plt
import seaborn as sns
from threadpoolctl import threadpool_limits
from sklearn.linear_model import LogisticRegressionCV
import base64
import io
import numpy as np
import scipy.stats as stats

from celery import shared_task
import pandas as pd
import decoupler as dc
import shap
import json
from .utils import (
    normalize1,
    loadSharedData,
    integrateCliData,
    integrateExData,
    getTopGeneCSV,
    clusteringPostProcess,
    go_it,
    delete_inactive_accounts,
    ICA,
    expression_clinic_split,
    normalize_string,
    map_normalized_keys,
    drop_parallel,
    MLparamSetting,
    shared_ml_dml_cf
)
from django.db import transaction
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import umap.umap_ as umap
from .models import MetaFileColumn, DatasetManager
import hdbscan
from sklearn.cluster import KMeans
import plotly.graph_objects as go
from sicaomics.annotate import toppfun

from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.metrics import roc_auc_score,mean_squared_error, roc_curve
import xgboost as xgb
import lightgbm as lgb


@shared_task
def cleanup_inactive_users():
    delete_inactive_accounts()



@shared_task(time_limit=180, soft_time_limit=150)
def runIntegrate(request, cID, usr):
    files = loadSharedData(request, cID)
    temp0 = integrateCliData(request, cID)

    if temp0.shape == (0, 0):
        raise Exception("Can't find meta file")
    if len(files) == 0:
        raise Exception("Can't find expression file")
    dfs1 = integrateExData(files, temp0)
    if dfs1 is None:
        raise Exception("No matched data for meta and omics")
    # combine Ex and clinic data
    temp = dfs1.set_index("ID_REF").join(
        temp0.set_index("ID_REF"), how="inner"
    )
    temp["obs"] = temp.index.tolist()
    # temp.dropna(axis=1, inplace=True)
    usr.setIntegrationData(temp)
    usr.metagenes = pd.DataFrame()
    usr.metageneCompose = pd.DataFrame()
    if usr.save() is False:
        raise Exception("Error for creating user records.")
    try:
        with transaction.atomic():
            new_file_columns = []
            MetaFileColumn.objects.filter(user=request.user, cID=cID).delete()
            for cn in temp0.columns:
                if cn == "LABEL":
                    label = "1"
                else:
                    label = "0"
                if np.issubdtype(temp0[cn].dtype, np.number):
                    num_flag = "1"
                else:
                    num_flag = "0"
                temp_meta = MetaFileColumn(
                    user=request.user,
                    cID=cID,
                    colName=cn,
                    label=label,
                    numeric=num_flag,
                )
                if temp_meta is None:
                    raise Exception("MetaFileColumn create Failed.")
                else:
                    new_file_columns.append(temp_meta)
            MetaFileColumn.objects.bulk_create(new_file_columns)
    except Exception as e:
        raise Exception(f"Error for registering to DataBase. {str(e)}")
    return

@shared_task(time_limit=180, soft_time_limit=150)
def runDgea(clusters, adata, targetLabel, n_genes):
    try:
        if clusters == "default":
            with plt.rc_context():
                figure1 = io.BytesIO()
                figure2 = io.BytesIO()
                with threadpool_limits(limits=NUMBER_CPU_LIMITS, user_api="blas"):
                    if len(set(adata.obs["batch1"])) > 1:
                        sc.tl.rank_genes_groups(
                            adata, groupby="batch1", method="t-test"
                        )
                        sc.tl.dendrogram(adata, groupby="batch1")
                        sc.pl.rank_genes_groups_dotplot(
                            adata, n_genes=int(n_genes), show=False, color_map="bwr"
                        )
                        plt.savefig(
                            figure1,
                            format="svg",
                            bbox_inches="tight",
                        )
                        figure1.seek(0)
                        figure1 = base64.b64encode(figure1.getvalue()).decode("utf-8")
                    else:
                        figure1 = ""
                    if len(set(adata.obs[targetLabel])) > 1:
                        adata.obs[targetLabel] = adata.obs[targetLabel].astype(str)
                        sc.tl.rank_genes_groups(
                            adata, groupby=targetLabel, method="t-test"
                        )
                        sc.tl.dendrogram(adata, groupby=targetLabel)
                        sc.pl.rank_genes_groups_dotplot(
                            adata, n_genes=int(n_genes), show=False, color_map="bwr"
                        )
                        plt.savefig(
                            figure2,
                            format="svg",
                            bbox_inches="tight",
                        )
                        figure2.seek(0)
                        figure2 = base64.b64encode(figure2.getvalue()).decode("utf-8")
                    else:
                        figure2 = ""
                plt.close('all')
                return [figure1, figure2]
        elif clusters == "fileName":
            return getTopGeneCSV(adata, "batch1", n_genes)
        elif clusters == "label":
            return getTopGeneCSV(adata, targetLabel, n_genes)
        elif clusters in ("LEIDEN", "HDBSCAN", "KMeans"):
            return getTopGeneCSV(adata, "cluster", n_genes)
    except:
        raise Exception("Error for running Dgea.")


@shared_task(time_limit=180, soft_time_limit=150)
def runClustering(cluster, adata, X2D, usr, param):
    if cluster == "LEIDEN":
        if param is None:
            param = 1
        try:
            param = float(param)
        except:
            raise Exception("Resolution should be a float")
        try:
            with threadpool_limits(limits=NUMBER_CPU_LIMITS, user_api="blas"):
                sc.pp.neighbors(adata, n_neighbors=40, n_pcs=40)
                sc.tl.leiden(adata, resolution=param)
        except Exception as e:
            raise Exception(f"{e}")
        Resp = clusteringPostProcess(X2D, adata, "leiden", usr)
        return Resp
    elif cluster == "HDBSCAN":
        if param is None:
            param = 20
        try:
            param = int(param)
        except:
            raise Exception("K should be positive integer.")
        if param < 5:
            raise Exception("minSize should be at least 5.")
        labels = hdbscan.HDBSCAN(min_cluster_size=int(param)).fit_predict(
            adata.obsm["X_pca"]
        )
        if min(labels) >= 0:
            labels = [str(i) for i in labels]
        else:
            labels = [str(i + 1) for i in labels]  # for outlier, will be assigned as -1

        adata.obs["hdbscan"] = labels
        adata.obs["hdbscan"] = adata.obs["hdbscan"].astype("category")
        Resp = clusteringPostProcess(X2D, adata, "hdbscan", usr)
        return Resp
    elif cluster == "Kmeans":
        try:
            param = int(param)
        except:
            raise Exception("K should be positive integer.")

        if param <= 1:
            raise Exception("K should be larger than 1.")
        km = KMeans(n_clusters=int(param), random_state=42, n_init="auto").fit(
            adata.obsm["X_pca"]
        )
        labels = [str(i) for i in km.labels_]
        adata.obs["kmeans"] = labels
        adata.obs["kmeans"] = adata.obs["kmeans"].astype("category")
        Resp = clusteringPostProcess(X2D, adata, "kmeans", usr)
        return Resp



@shared_task(time_limit=180, soft_time_limit=150)
def runFeRed(fr, usr):
    pca_temp = usr.getAnndata().obsm["X_pca"]
    if fr == "TSNE":
        tsne = TSNE(
            n_components=2,
            random_state=42,
            n_jobs=2,
            perplexity=min(30.0, pca_temp.shape[0] - 1),
        )
        X2D = tsne.fit_transform(pca_temp)
        usr.redMethod = "TSNE"
    elif fr == "UMAP":
        umap1 = umap.UMAP(
            n_components=2,
            random_state=42,
            n_neighbors=min(30, pca_temp.shape[0] // 2),
            metric="euclidean",
            init="spectral",
            n_epochs=None,
            learning_rate=1.0,
            n_jobs=2,
            spread=1.0,
            min_dist=0.1,
        )
        X2D = umap1.fit_transform(pca_temp)
        usr.redMethod = "UMAP"
    else:
        pca = PCA(n_components=2, random_state=42)
        X2D = pca.fit_transform(pca_temp)
        usr.redMethod = "PC"
    return X2D


@shared_task(time_limit=180, soft_time_limit=150)
def runICA(df, num, suffix='', flag1=None, flag2=None):
    df = df.copy()
    df.index.name = "ID_REF"
    # Drop the columns 'FileName', 'obs', and 'LABEL'

    columns_to_drop = ["FileName", "obs", "LABEL"]
    df.drop(columns=[col for col in columns_to_drop if col in df.columns], inplace=True)
    df1 = df.loc[:, ~df.columns.str.startswith("c_")]
    # cell_type is compulsory
    if 'cell_type' in df1.columns:
        df1 = pd.get_dummies(df1, columns=['cell_type'], drop_first=False)

    try:
        if flag1 is None:
            metageneCompose, metagenes, _, sICA = ICA(df1, None, num, '')
            metagenes2 = None
        else:
            df11 = df1[flag1]
            df12 = df1[flag2]
            metageneCompose, metagenes, metagenes2, sICA = ICA(df11, df12, num, suffix)
    except Exception as e:
        raise e
    return metageneCompose, metagenes, metagenes2,sICA
    
@shared_task(time_limit=180, soft_time_limit=150)
def runPCA(df, num):
    pca = PCA(n_components=num)
    df_pca = pca.fit_transform(df)
    metagenes = pd.DataFrame(
        df_pca, 
        columns = ['PCA'+str(i) for i in range(num)],
        index=df.index
    )
    return metagenes


@shared_task(time_limit=180, soft_time_limit=150)
def runTopFun(metageneCompose, thre=3):
    Tannot = toppfun.ToppFunAnalysis(
        data=metageneCompose,
        threshold=thre,
        method="std",
        tail="heaviest",
        convert_ids=True,
    )
    
    return Tannot.top_genes_[['inputs']]

@shared_task(time_limit=180, soft_time_limit=150)
def calculate_imf_scores(expression_data):
    gene_sets = {
        'I': {  # Immune cells
            'CD8_T': ['CD8A', 'CD8B', 'CD3D', 'CD3E', 'CD3G'],
            'CD4_T': ['CD4', 'CD3D', 'CD3E', 'CD3G'],
            'B_cells': ['CD19', 'CD79A', 'CD79B', 'MS4A1'],
            'NK': ['NCAM1', 'KLRB1', 'KLRD1', 'NKG7']
        },
        'M': {  # Myeloid cells
            'Monocytes': ['CD14', 'CD163', 'CSF1R'],
            'Macrophages': ['CD68', 'CD163', 'CSF1R'],
            'Neutrophils': ['CEACAM8', 'FCGR3B', 'CSF3R'],
            'DC': ['CD1C', 'CLEC4C', 'NRP1']
        },
        'F': {  # Fibroblasts
            'Fibroblasts': ['COL1A1', 'COL1A2', 'FAP', 'PDGFRB'],
            'Myofibroblasts': ['ACTA2', 'TAGLN', 'MYL9']
        }
    }
    # Initialize results DataFrame
    imf_scores = pd.DataFrame(index=expression_data.index, columns=['I', 'M', 'F'])

    # Calculate scores for each component
    for component, cell_types in gene_sets.items():
        component_scores = []
        
        for cell_type, genes in cell_types.items():
            # Find which genes from our signature are present in the data
            available_genes = [gene for gene in genes if gene in expression_data.columns]
            
            if len(available_genes)>0:
                # Calculate mean expression for this cell type
                cell_score = expression_data[available_genes].mean(axis=1)
                component_scores.append(cell_score)
            else:
                raise Exception(f"Warning: No genes found for {cell_type} in {component} score calculation")
        
        if len(component_scores) > 0:
            # Average the scores for all cell types in this component
            imf_scores[component] = np.mean(component_scores, axis=0)
        else:
            raise Exception(f"Warning: No scores could be calculated for {component}")
    
    # Normalize scores to 0-1 range
    for col in ['I', 'M', 'F']:
        min_val, max_val = imf_scores[col].min(), imf_scores[col].max()
        if max_val > min_val:  # Avoid division by zero
            imf_scores[col] = (imf_scores[col] - min_val) / (max_val - min_val)
        else:
            raise Exception(f"imf min==max.")
    
    def fig_add_imf(imf_scores):
        pfig, axes = plt.subplots(1, 3, figsize=(12, 4))

        for i, component in enumerate(['I', 'M', 'F']):
            sns.histplot(imf_scores[component], kde=True, ax=axes[i], color=["#1f77b4", "#ff7f0e", "#2ca02c"][i])
            axes[i].set_title(f'{component} Score Distribution')
            axes[i].set_xlabel('Score')

        plt.tight_layout()
        encoded_img = base64.b64encode(fig_to_image(pfig, "svg")).decode("utf-8")
        plt.close(pfig)
        return encoded_img
        
    def fig_to_image(fig, fmt="svg"):
        from io import BytesIO
        buf = BytesIO()
        fig.savefig(buf, format=fmt)
        buf.seek(0)
        return buf.read()
        
    fig = fig_add_imf(imf_scores)
    return imf_scores,fig

import decoupler as dc
@shared_task(time_limit=60*10, soft_time_limit=580)
def run_ora_task(db, adata1,de,msigdb1):
    try:
        dc.run_ora(
                mat=adata1,
                net=msigdb1,
                source='geneset',
                target='genesymbol',
                verbose=True,
                use_raw=False
        )
    except Exception as e:
        raise Exception(f"Error during ORA computation: {e}")
    acts = dc.get_acts(adata1, obsm_key='ora_estimate')
    acts_v = acts.X.ravel()
    max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
    acts.X[~np.isfinite(acts.X)] = max_e
    dff1 = dc.rank_sources_groups(acts, groupby=de, reference='rest', method='t-test_overestim_var')[['group','reference','names','pvals','pvals_adj']]
    return dff1,acts


@shared_task(time_limit=60*3, soft_time_limit=160)
def run_oraPlot(acts, source_markers, group):
    fea = list(set().union(*source_markers.values()))
    if group=='LABEL':
        group='batch2'
    sc.set_figure_params(dpi=100)
    sc.settings.verbosity = 0
    # Iterate over genes and plot
    with threadpool_limits(limits=NUMBER_CPU_LIMITS, user_api="blas"):
        figure1 = io.BytesIO()
        sc.pl.matrixplot(acts, fea, group, dendrogram=True, standard_scale='var', colorbar_title='Z-scaled scores', cmap='RdBu_r')
        #sns.violinplot(data=acts, x="group", y="pvals_adj", inner="point", palette="coolwarm")
        plt.savefig(figure1, format="svg", bbox_inches="tight")
        plt.close()
    image_data = base64.b64encode(figure1.getvalue()).decode("utf-8")
    return image_data

@shared_task(time_limit=60*4, soft_time_limit=230)
def run_Rf(typeRf, param, X, X_train, y_train, X_test, y_test, testRatio):
    dataset_manager = DatasetManager(X, X_train, X_test)
    if typeRf=='classification':
        rf = RandomForestClassifier(**param)
        rf.fit(X_train, y_train)
        y_prob = rf.predict_proba(X_test)[:, 1]  # Positive class probabilities
        auc_score = roc_auc_score(y_test, y_prob)
        returnDS = dataset_manager.get_dataset(SHAP_PLOT_DATASET)
        return auc_score, (rf, returnDS)     
    else:
        rf = RandomForestRegressor(**param)
        rf.fit(X_train, y_train)
        y_pred = rf.predict(X_test)
        mse = mean_squared_error(y_test, y_pred)
        returnDS = dataset_manager.get_dataset(SHAP_PLOT_DATASET)
        return mse, (rf, returnDS)
        
@shared_task(time_limit=60*4, soft_time_limit=230)
def run_XGBoost(typeRf, param, X, X_train, y_train, X_test, y_test, testRatio):
    dataset_manager = DatasetManager(X, X_train, X_test)
    if typeRf == 'classification':
        model = xgb.XGBClassifier(**param)
        model.fit(X_train, y_train)
        y_prob = model.predict_proba(X_test)[:, 1]  # Positive class probabilities
        auc_score = roc_auc_score(y_test, y_prob)
        returnDS = dataset_manager.get_dataset(SHAP_PLOT_DATASET)
        return auc_score, (model, returnDS)
    else:
        model = xgb.XGBRegressor(**param)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        mse = mean_squared_error(y_test, y_pred)
        returnDS = dataset_manager.get_dataset(SHAP_PLOT_DATASET)
        return mse, (model, returnDS)

@shared_task(time_limit=60*4, soft_time_limit=230)
def run_LightGBM(typeRf, param, X, X_train, y_train, X_test, y_test, testRatio):
    dataset_manager = DatasetManager(X, X_train, X_test)
    if typeRf == 'classification':
        model = lgb.LGBMClassifier(**param)
        model.fit(X_train, y_train)
        y_prob = model.predict_proba(X_test)[:, 1]  # Positive class probabilities
        auc_score = roc_auc_score(y_test, y_prob)
        returnDS = dataset_manager.get_dataset(SHAP_PLOT_DATASET)
        return auc_score, (model, returnDS)
    else:
        model = lgb.LGBMRegressor(**param)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        mse = mean_squared_error(y_test, y_pred)
        returnDS = dataset_manager.get_dataset(SHAP_PLOT_DATASET)
        return mse, (model, returnDS)
        
@shared_task(time_limit=60*4, soft_time_limit=230)
def run_GB(typeRf, param, X, X_train, y_train, X_test, y_test, testRatio):
    dataset_manager = DatasetManager(X, X_train, X_test)
    if typeRf == 'classification':
        model = GradientBoostingClassifier(**param)
        model.fit(X_train, y_train)
        y_prob = model.predict_proba(X_test)[:, 1]  # Positive class probabilities
        auc_score = roc_auc_score(y_test, y_prob)
        returnDS = dataset_manager.get_dataset(SHAP_PLOT_DATASET)
        return auc_score, (model, returnDS)
    else:
        model = GradientBoostingRegressor(**param)
        model.fit(X_train, y_train)

        # Predict and calculate MSE
        y_pred = model.predict(X_test)
        mse = mean_squared_error(y_test, y_pred)

        returnDS = dataset_manager.get_dataset(SHAP_PLOT_DATASET)
        return mse, (model, returnDS)
        
@shared_task(time_limit=60*4, soft_time_limit=230)
def run_shap(model, X, flip='no'):
    explainer = shap.TreeExplainer(model, X)
    shap_values = explainer(X, check_additivity=False)
    if len(shap_values.shape)!=2:
        shap_values = shap_values[:,:,1]
    feature_order = None
    if flip=='yes':
        mean_shap_importance = np.abs(shap_values.values).mean(axis=0)
        feature_order = np.argsort(mean_shap_importance)
    with threadpool_limits(limits=NUMBER_CPU_LIMITS, user_api="blas"):
        figure1 = io.BytesIO()
        if flip=='no':
            shap.plots.beeswarm(shap_values, max_display=40)
        else:
            shap.plots.beeswarm(shap_values, max_display=40, show=False, order=feature_order)
        plt.savefig(figure1, format="svg", bbox_inches="tight")
        plt.close()
    image_data = base64.b64encode(figure1.getvalue()).decode("utf-8")
    return image_data, shap_values

@shared_task(time_limit=80, soft_time_limit=70)    
def run_MLpreprocess(label, cla, testRatio, logY, usr, mlMethod, mlType,  ml_params=''):
    if label=='LABEL':
        label='batch2'
    testRatio = 0.1 if testRatio == '1' else 0.2
    param = MLparamSetting(mlMethod, ml_params)
    if label is None:
        raise Exception("incorrect message request")
    X, y, T = shared_ml_dml_cf(usr, logY, mlType, label, cla)
    return mlType, param, X, y, testRatio, T
    
@shared_task(time_limit=600, soft_time_limit=580)
def run_DMLreport(model1, X1, model2, X2):
    def compute_shap_importance(model, X):
        explainer = shap.TreeExplainer(model, X)
        shap_values = explainer(X)
        if len(shap_values.shape)!=2:
            shap_values = shap_values[:,:,1]
        shap_df = pd.DataFrame(shap_values.values, columns=X.columns)
        return shap_values, shap_df.abs().mean(axis=0)
        
    shap_values1, shap_importance1 = compute_shap_importance(model1, X1)
    shap_values2, shap_importance2 = compute_shap_importance(model2, X2)
    
    shap_importance3 = shap_importance1 * shap_importance2
    shap_importance_sorted = shap_importance3.sort_values(ascending=False)
    
    shap_min = shap_importance_sorted.min()
    shap_max = shap_importance_sorted.max()
    
    with threadpool_limits(limits=NUMBER_CPU_LIMITS, user_api="blas"):
        fig, ax = plt.subplots()
        shap_importance_sorted.plot(kind='barh', color='steelblue', ax=ax)
        ax.set_xlim(shap_min, shap_max)  # Ensure consistent x-axis across all plots
        ax.invert_yaxis()  # Highest values on top
        ax.set_xlabel("Mean Absolute SHAP Value")
        ax.set_ylabel("Features")
        ax.set_title(f"SHAP importance on cohorts")
        ax.set_aspect('auto')
        plt.tight_layout()
        figure1 = io.BytesIO()
        plt.savefig(figure1, format="svg", bbox_inches="tight")
        plt.close()
    image_data = base64.b64encode(figure1.getvalue()).decode("utf-8")
    return image_data, (shap_values1, shap_values2), shap_importance_sorted.head(6).index.tolist()
    
@shared_task(time_limit=40, soft_time_limit=30)
def run_DMLreportDis(metagenes, shap_importance_sorted1, shap_importance_sorted2, ica_cohort):
    num_plots = len(metagenes)
    rows = (num_plots + 1) // 2
    fig, axes = plt.subplots(rows, 2, figsize=(10, rows * 3))
    axes = axes.flatten()  # Flatten for easier indexing

    # Determine dynamic SHAP limits
    max_shap = max(shap_importance_sorted1.max(), shap_importance_sorted2.max())
    shap_limit = max_shap * 1.1  # Add margin for visibility

    with threadpool_limits(limits=NUMBER_CPU_LIMITS, user_api="blas"):
        for k, gene in enumerate(metagenes):
            ax = axes[k]

            # Extract SHAP values safely
            val1 = shap_importance_sorted1.get(gene, 0)
            val2 = shap_importance_sorted2.get(gene, 0)

            # Plot SHAP values (downward/upward)
            ax.bar('val1', val1, color='red', alpha=0.7, label=ica_cohort[1])
            ax.bar('val2', val2, color='green', alpha=0.7, label=ica_cohort[2])

            # Dynamic axis limits
            ax.set_ylim(-shap_limit, shap_limit)

            # Set x-axis and title
            ax.set_xticks(['val1', 'val2'])
            ax.set_xticklabels([gene, gene])
            ax.set_title(f"SHAP on coefficient: {gene}")

            if (k % 2) == 1:  # Second column (0-based indexing)
                ax.set_ylabel('')
            else:
                ax.set_ylabel("SHAP Importance Value")
            ax.legend()

    # Handle cases where there are fewer plots than grid slots
    for ax in axes[num_plots:]:
        fig.delaxes(ax)
    plt.tight_layout(pad=2.0)
    # Save figure to memory
    figure1 = io.BytesIO()
    plt.savefig(figure1, format="svg", bbox_inches="tight")
    plt.close()

    # Encode image and return JSON response
    image_data = base64.b64encode(figure1.getvalue()).decode("utf-8")
    return image_data
    
@shared_task(time_limit=40, soft_time_limit=30)
def run_CFreport1(mlType, model_t, model_y, X_test, T_test, y_test):
    def plot_model_performance(yType, model_t, model_y, X_test, T_test, Y_test):
        figure = io.BytesIO()
        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1)
        T_pred = model_t.predict_proba(X_test)[:, 1]
        plot_roc_curve(T_test, T_pred, "Treatment Model")
        plt.subplot(1, 2, 2)
        if yType == 'regression':
            Y_pred = model_y.predict(X_test)
            mse = mean_squared_error(Y_test, Y_pred)
            plot_pred_vs_actual(Y_test, Y_pred, "Outcome Model", mse)
        else:
            Y_pred = model_y.predict_proba(X_test)[:, 1]
            plot_roc_curve(Y_test, Y_pred, "Outcome Model")
        plt.tight_layout()
        plt.savefig(figure, format="svg", bbox_inches="tight")
        plt.close()
        return base64.b64encode(figure.getvalue()).decode("utf-8")

    def plot_roc_curve(y_true, y_pred, title):
        fpr, tpr, _ = roc_curve(y_true, y_pred)
        auc = roc_auc_score(y_true, y_pred)
        plt.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(f'{title} ROC Curve')
        plt.legend()

    def plot_pred_vs_actual(y_true, y_pred, title, mse=None):
        plt.scatter(y_true, y_pred, alpha=0.5)
        plt.plot([min(y_true), max(y_true)], [min(y_true), max(y_true)], 'r--')
        plt.xlabel("Actual Y")
        plt.ylabel("Predicted Y")
        plt.title(f"{title}: Predicted vs Actual")
        if mse is not None:
            plt.text(0.1, 0.9, f'MSE: {mse:.4f}', ha='center', va='center', transform=plt.gca().transAxes, fontsize=12, color='blue')
        plt.grid()
    
    image = plot_model_performance(mlType, model_t, model_y, X_test, T_test, y_test)
    return image
    
@shared_task(time_limit=40, soft_time_limit=30)
def run_CFreport2(causal_forest, X_test):
    treatment_effects = causal_forest.effect(X_test)
    sns.histplot(treatment_effects.flatten(), bins=30, kde=True)
    figure = io.BytesIO()
    plt.axvline(0, color='red', linestyle='--')  # Mark neutral effect
    plt.xlabel(r"Estimated Treatment Effect ($\hat{\tau}(X)$)")
    plt.ylabel("Count")
    plt.title("Distribution of Estimated Treatment Effects")
    plt.savefig(figure, format="svg", bbox_inches="tight")
    plt.close()
    return base64.b64encode(figure.getvalue()).decode("utf-8")
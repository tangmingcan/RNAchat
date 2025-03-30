from django.urls import path, re_path
from . import views

urlpatterns = [
    # path("", views.index, name="index"),
    path("accounts/rest/login/", views.restLogin, name="restLogin"),
    path("geneExpression/", views.opExpression, name="op_expression"),
    path("meta/", views.opMeta, name="op_meta"),
    path("", views.tab, name="tab"),
    path("meta/columns/", views.meta_columns, name="meta_columns"),
    path("eda/integrate/", views.edaIntegrate, name="edaIntegrate"),
    path(
        r"meta/<slug:colName>/",
        views.meta_column_values,
        name="meta_columns_values",
    ),
    path("user/", views.checkUser, name="checkUser"),
    path("ica/metagenes/", views.downloadICA, name="downloadICA"),
    path("compType/ICA/report/", views.ICAreport, name="ICAreport"),
    path("compType/ICA/ica_cohort/", views.compICAcohort, name="compICAcohort"),
    path("compType/ICA/ica_all/", views.compICAall, name="compICAall"),
    path("ml/featureImportance/data/", views.feaImpData, name="feaImpData"),
    path("ml/featureImportance/", views.feaImp, name="feaImp"),
    path("ml/featureDependence/", views.ShapInteraction, name="ShapInteraction"),
    path("ml/<slug:mlType>/<slug:mlMethod>/", views.goML, name="runML"),
    path("dml/featureImportance/<int:number>/", views.DMLfeaImp, name="DMLfeaImp"),
    path("dml/report/", views.DMLreport, name="DMLreport"),
    path("dml/report/distribution/",views.DMLreportDis, name="DMLreportDis"),
    path("dml/<slug:mlType>/<slug:mlMethod>/<int:number>/", views.goDML, name="runDML"),
    path("cf/report/TYmodels/", views.CFreport1, name="CFreport1"),
    path("cf/report/causalForests/", views.CFreport2, name="CFreport2"),
    path("cf/featureImportance/", views.CFfeaImp, name="CFfeaImp"),
    path("dataVisual1/", views.dataVisual1, name="dataVisual1"),
    path("dataVisual1/matrixplot/", views.matrixplot, name="matrixplot"),
    path("dataVisual2/metagenes/", views.dataVisual2, name="dataVisual2")
]

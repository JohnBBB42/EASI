# easi_analysis/urls.py

from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('cv/', views.cv_analysis_view, name='cv_analysis'),
    path('dpv/', views.dpv_analysis_view, name='dpv_analysis'),
    path('eis/', views.eis_analysis_view, name='eis_analysis'),
]

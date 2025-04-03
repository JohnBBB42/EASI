# easi_analysis/urls.py
from django.conf import settings
from django.conf.urls.static import static
from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('cv/', views.cv_analysis_view, name='cv_analysis'),
    path('dpv/', views.dpv_analysis_view, name='dpv_analysis'),
    path('eis/', views.eis_analysis_view, name='eis_analysis'),
]

# ✅ Add this block to serve static files during dev
if settings.DEBUG:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

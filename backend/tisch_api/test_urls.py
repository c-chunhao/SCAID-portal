from django.urls import include, path

urlpatterns = [path('api/', include('tisch_api.urls'))]

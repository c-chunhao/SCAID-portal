class CSPMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        response = self.get_response(request)
        response["Content-Security-Policy"] = (
            "img-src 'self' data: blob: http://localhost:8000/ "
            "default-src 'self';"
        )
        return response
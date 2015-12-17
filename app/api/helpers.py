from flask.globals import _app_ctx_stack, _request_ctx_stack
from werkzeug.urls import url_parse
from werkzeug.exceptions import NotFound


def match_url(url, method=None):
    """Get endpoint and arguments from a url.

    SEE http://stackoverflow.com/questions/19631335/reverting-a-url-in-flask

    Returns:
        endpoint name (str?)
        arguments (dict).
    """
    appctx = _app_ctx_stack.top
    reqctx = _request_ctx_stack.top
    if appctx is None:
        raise RuntimeError('Attempted to match a URL without the '
                           'application context being pushed. This has to be '
                           'executed when application context is available.')

    if reqctx is not None:
        url_adapter = reqctx.url_adapter
    else:
        url_adapter = appctx.url_adapter
        if url_adapter is None:
            raise RuntimeError('Application was not able to create a URL '
                               'adapter for request independent URL matching. '
                               'You might be able to fix this by setting '
                               'the SERVER_NAME config variable.')
    parsed_url = url_parse(url)
    if parsed_url.netloc is not '' and \
            parsed_url.netloc != url_adapter.server_name:
        raise NotFound()
    return url_adapter.match(parsed_url.path, method)


def args_from_url(url, endpoint):
    """Get argument dictionary from url, if endpoint is correct."""
    r = match_url(url, 'GET')
    if r[0] != endpoint:
        raise NotFound()
    return r[1]

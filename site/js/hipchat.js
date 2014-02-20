(function() {

    // Creates an iframe with an embedded HipChat conversation window.
    //
    // Options:
    //     url - The url to the room to embed; required
    //     el - The container in which to insert the HipChat panel; required
    //     timezone - The timezone to use in the embedded room; required
    //     width - The width of the iframe; defaults to 100%
    //     height - The height of the iframe; defaults to 400px

    var parametize = function(params) {
        var key, toks = [];
        for (key in params) {
            toks.push(key + '=' + params[key]);
        }
        return toks.join('&');
    };

    return this.HipChat = function(options) {
        if (options && options.url && options.el && options.timezone) {
            var el = document.querySelector(options.el);
            if (!el) return;
            var params = {
                anonymous: 1,
                timezone: options.timezone,
                minimal: 1
            };
            var url = options.url + (options.url.indexOf('?') > 0 ? '&' : '?') +
                parametize(params);
            if (url.indexOf('https://') !== 0) {
                url = 'https://' + url;
            }
            var w = options.width || '100%';
            var h = options.height || 600;
            el.innerHTML = '<iframe src="' + url + '" frameborder="' + 0 +
                '" width="' + w + '" height="' + h + '"></iframe>';
        }
    };
})();

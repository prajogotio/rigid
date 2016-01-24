function include(src, state) {
	var script = document.createElement('script');
	script.type = 'text/javascript';
	script.src = src;
	script.onload = function() {
		state.includes++;
		if (state.includes == state.NUM_OF_INCLUDES) {
			state.run();
		}
	}
	document.getElementsByTagName('head')[0].appendChild(script);
}

function includeAllScripts() {
	var state = {
		src : ['Vec2.js', 'Body.js', 'Sphere.js', 'Manifold.js'],
		includes : 0,
		run : function() {
			testPlatform();
		}
	}
	state.NUM_OF_INCLUDES = state.src.length;

	for (var i = 0; i < state.NUM_OF_INCLUDES; ++i) {
		include(state.src[i], state);
	}
}


addEventListener("DOMContentLoaded", function() {
	includeAllScripts();
});
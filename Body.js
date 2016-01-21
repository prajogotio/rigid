// Must be convex
function Body(cm, arrayOfVertices, 
			  mass, e, momentOfInertia,
			  velocity, acceleration,
			  angularVelocity) {
	this.v = arrayOfVertices;	// set of vertices
	this.n = [];				// normals of edges defined by v[i] and v[i+1]
	for (var i = 0; i < this.v.length; ++i) {
		var j = i+1;
		if (j == this.v.length) j = 0;
		this.n.push(getNormal(this.v[i], this.v[j]));
	}

	this.cm = cm; // center of mass, Vec2
	this.velocity = velocity || new Vec2(0, 0); // velocity, Vec2
	this.acceleration = acceleration || new Vec2(0, 0); // acceleration

	this.mass = mass;	// scalar
	this.inverseMass = (mass == 0 ? 0 : 1/mass);
	this.e = e; 		// scalar resitution of elasticity

	this.momentOfInertia = momentOfInertia;	// scalar
	this.inverseMomentOfInertia = (momentOfInertia == 0 ? 0 : 1/momentOfInertia);
	this.angularVelocity = angularVelocity || 0; // scalar radian
	this.angularPosition = 0;	// radian

}

Body.prototype.getAbsolutePosition = function(v) {
	return v.copy().rotate(this.angularPosition/Math.PI*180).plus(this.cm);
}

Body.prototype.getAbsoluteDirection = function(n) {
	return n.copy().rotate(this.angularPosition/Math.PI*180);
}

Body.prototype.axisOfLeastSeparationWith = function(b) {
	// b is another Body
	var MINIMUM_LENGTH = -1e301;
	var MAXIMUM_LENGTH = 1e301;
	var chosen = -1;
	var minLength = MAXIMUM_LENGTH;
	for (var i = 0; i < this.v.length; ++i) {
		var curMaxLength = MINIMUM_LENGTH;
		for (var j = 0; j < b.v.length; ++j) {
			var cur = this.getAbsoluteDirection(this.n[i]).flip().dot(b.getAbsolutePosition(b.v[j]).minus(this.getAbsolutePosition(this.v[i])));
			curMaxLength = Math.max(curMaxLength, cur);
		}
		if (curMaxLength < minLength) {
			minLength = curMaxLength;
			chosen = i;
		}
	}
	return {
		axis: this.getAbsoluteDirection(this.n[chosen]),
		depth: minLength,
	};
}

Body.prototype.getBestContactEdge = function(n) {
	var MINIMUM_DISTANCE = -1e301;
	var chosen = -1;
	var furthest = MINIMUM_DISTANCE;
	for (var i = 0; i < this.v.length; ++i) {
		var cur = this.getAbsolutePosition(this.v[i]).dot(n);
		if (furthest < cur) {
			furthest = cur;
			chosen = i;
		}
	}
	var prev = chosen-1;
	var next = (chosen+1)%this.v.length;
	if (prev < 0) prev += this.v.length;
	var pval = Math.abs(n.dot(this.getAbsolutePosition(this.v[prev]).minus(this.getAbsolutePosition(this.v[chosen]))));
	var nval = Math.abs(n.dot(this.getAbsolutePosition(this.v[chosen]).minus(this.getAbsolutePosition(this.v[next]))));
	// edge [vec2, vec2] preserves relative direction as in the convex polygon
	// so that the function getNormal returns the correct normal
	return (pval < nval ? [this.getAbsolutePosition(this.v[prev]), this.getAbsolutePosition(this.v[chosen])]
						: [this.getAbsolutePosition(this.v[chosen]), this.getAbsolutePosition(this.v[next])]);
}

Body.prototype.collidesWith = function(b) {
	// separating axis theory
	for (var i = 0; i < this.v.length; ++i) {
		if (!checkProjectionsOverlap(this, b, this.getAbsoluteDirection(this.n[i]))) {
			return false;
		}
	}
	for (var i = 0; i < b.v.length; ++i) {
		if (!checkProjectionsOverlap(this, b, b.getAbsoluteDirection(b.n[i]))) {
			return false;
		}
	}
	return true;
}

Body.prototype.integrate = function(dt) {
	// translational
	this.velocity = this.velocity.plus(this.acceleration.times(dt));
	this.cm = this.cm.plus(this.velocity.times(dt));

	// rotational
	this.angularPosition += this.angularVelocity;
}

Body.prototype.applyImpulse = function(j, n) {
	// j is scalar amount of impulse
	// n is the direction of impulse

	// translational
	this.velocity = this.velocity.plus(n.times(j*this.inverseMass));
}


function getContactPoints(ref, inc) {
	// ref reference edge,
	// inc incident edge

	var contacts = [inc[0].copy(), inc[1].copy()];
	var n = getNormal(ref[0], ref[1]);
	
	// first clipping
	contacts = clip(ref[0], ref[1].minus(ref[0]), contacts, n, true);
	// second clipping
	contacts = clip(ref[1], ref[0].minus(ref[1]), contacts, n, true);

	// third clipping: without replacement
	contacts = clip(ref[0], n.flip(), contacts, n, false);
	return contacts;
}

function clip(A, dir, points, n, withReplacement) {
	// clipping with replacement if available
	// points must contain 2 points (an edge)
	var ret = [];
	var outside = null;	// the point that is clipped
	for (var i = 0; i < points.length; ++i) {
		if (halfPlaneContains(A, dir, points[i])) {
			ret.push(points[i].copy());
		} else {
			outside = points[i].copy();
		}
	}

	if (ret.length == 0 || ret.length == points.length) return ret;

	if (withReplacement){
		// invariance: points must contain 2 points, ret must contain 1 point
		var res = rayWithSegmentIntersection(A, n, ret[0], outside);
		ret.push(res.intersect);
	}
	return ret;
}

function checkProjectionsOverlap(A, B, n) {
	var MAXIMUM_VALUE = 1e301, MINIMUM_VALUE = -1e301;
	var Amin = MAXIMUM_VALUE,
		Amax = MINIMUM_VALUE,
		Bmin = MAXIMUM_VALUE,
		Bmax = MINIMUM_VALUE;
	// n should have been normalized
	n.normalize();
	for (var i = 0; i < A.v.length; ++i) {
		var cur = A.getAbsolutePosition(A.v[i]).dot(n);
		Amin = Math.min(Amin, cur); 
		Amax = Math.max(Amax, cur);
	}
	for (var i = 0; i < B.v.length; ++i) {
		var cur = B.getAbsolutePosition(B.v[i]).dot(n);
		Bmin = Math.min(Bmin, cur);
		Bmax = Math.max(Bmax, cur);
	}
	return !(Amax <= Bmin || Bmax <= Amin);
}

function resolveCollision(A, B) {
	// polygon-polygon resolver for point B
	// A is the reference


	// choose a normal from A
	// n must point out of A
	var res = A.axisOfLeastSeparationWith(B);
	var n = res.axis;
	var depth = res.depth;

	// getting contact of incident B
	var ref = A.getBestContactEdge(n);
	var inc = B.getBestContactEdge(n.flip());
	var contacts = getContactPoints(ref, inc);


	// relative velocity of B as seen from A
	var vrel = B.velocity.minus(A.velocity);
	var proj = vrel.dot(n);
	if (proj > 1e-12) {
		// objects moving away from each other
		return;
	}
	if (proj > -1e-12) {
		positionalCorrection(A, B, n, depth);
		return;
	}
	// use the one with lower elasticity
	var e = Math.min(A.e, B.e);
	var j = -(1+e)*proj/(A.inverseMass + B.inverseMass);
	A.applyImpulse(-j, n);
	B.applyImpulse(j, n);
	positionalCorrection(A, B, n, depth);
}

function positionalCorrection(A, B, n, depth) {
	// A is the reference body
	// n points out of A
	var CORRECTION_RATE = 0.2;
	var correction = depth/(A.inverseMass + B.inverseMass);
	A.cm = A.cm.plus(n.flip().times(correction*A.inverseMass));
	B.cm = B.cm.plus(n.times(correction*B.inverseMass));
	
}
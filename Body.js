// Polygon body. Must be convex
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
	this.torque = 0;
	this.angularVelocity = angularVelocity || 0; // scalar radian
	this.angularPosition = 0;	// radian

	this.staticFriction = 0.3;
	this.dynamicFriction = 0.22;

	this.ANGULAR_DECAY = 0.98;

	this.cachedV = this.v.slice(0);
	this.cachedN = this.n.slice(0);
	this.boundingRectA = new Vec2(0,0);
	this.boundingRectB = new Vec2(0,0);

	this.recomputeModel();
}

Body.prototype.recomputeModel = function() {
	this.boundingRectA.x = this.boundingRectA.y = 1e301;
	this.boundingRectB.x = this.boundingRectB.y = -1e301;
	for (var i = 0; i < this.v.length; ++i) {
		this.cachedV[i] = this.v[i].copy().rotate(this.angularPosition/Math.PI*180).plus(this.cm);
		this.cachedN[i] = this.n[i].copy().rotate(this.angularPosition/Math.PI*180);

		this.boundingRectA.x = Math.min(this.boundingRectA.x, this.cachedV[i].x);
		this.boundingRectA.y = Math.min(this.boundingRectA.y, this.cachedV[i].y);
		this.boundingRectB.x = Math.max(this.boundingRectB.x, this.cachedV[i].x);
		this.boundingRectB.y = Math.max(this.boundingRectB.y, this.cachedV[i].y);
	}
}

Body.prototype.getAbsolutePosition = function(i) {
	return this.cachedV[i];
}

Body.prototype.getAbsoluteDirection = function(i) {
	return this.cachedN[i];
}

Body.prototype.axisOfLeastSeparationWith = function(b) {
	if (b instanceof Sphere) {
		return b.axisOfLeastSeparationWith(this);
	}

	// b is another Polygon
	var MINIMUM_LENGTH = -1e301;
	var MAXIMUM_LENGTH = 1e301;
	var chosen = -1;
	var minLength = MAXIMUM_LENGTH;
	for (var i = 0; i < this.v.length; ++i) {
		var curMaxLength = MINIMUM_LENGTH;
		for (var j = 0; j < b.v.length; ++j) {
			var cur = this.getAbsoluteDirection(i).flip().dot(b.getAbsolutePosition(j).minus(this.getAbsolutePosition(i)));
			curMaxLength = Math.max(curMaxLength, cur);
		}
		if (curMaxLength < minLength) {
			minLength = curMaxLength;
			chosen = i;
		}
	}

	// support point way
	var axis = this.getAbsoluteDirection(chosen);


	return {
		axis: axis,
		depth: minLength,
	};
}

Body.prototype.getBestContactEdge = function(n) {
	var MINIMUM_DISTANCE = -1e301;
	var chosen = 0;
	var furthest = MINIMUM_DISTANCE;
	for (var i = 0; i < this.v.length; ++i) {
		var cur = this.getAbsolutePosition(i).dot(n);
		if (furthest < cur) {
			furthest = cur;
			chosen = i;
		}
	}
	var prev = chosen-1;
	var next = (chosen+1)%this.v.length;
	if (prev < 0) prev += this.v.length;
	var pval = Math.abs(n.dot(this.getAbsolutePosition(prev).minus(this.getAbsolutePosition(chosen))));
	var nval = Math.abs(n.dot(this.getAbsolutePosition(chosen).minus(this.getAbsolutePosition(next))));
	// edge [vec2, vec2] preserves relative direction as in the convex polygon
	// so that the function getNormal returns the correct normal
	return (pval < nval ? [this.getAbsolutePosition(prev), this.getAbsolutePosition(chosen)]
						: [this.getAbsolutePosition(chosen), this.getAbsolutePosition(next)]);
}

Body.prototype.checkBoundingRect = function(b) {
	return !(this.boundingRectB.x < b.boundingRectA.x || b.boundingRectB.x < this.boundingRectA.x
		|| this.boundingRectB.y < b.boundingRectA.y || b.boundingRectB.y < this.boundingRectA.y);
}

Body.prototype.collidesWith = function(b) {

	if (b instanceof Sphere) {
		return b.collidesWith(this);
	}
	if (!this.checkBoundingRect(b)) return false;
	// separating axis theory
	for (var i = 0; i < this.v.length; ++i) {
		if (!checkProjectionsOverlap(this, b, this.getAbsoluteDirection(i))) {
			return false;
		}
	}
	for (var i = 0; i < b.v.length; ++i) {
		if (!checkProjectionsOverlap(this, b, b.getAbsoluteDirection(i))) {
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
	//this.angularVelocity += this.torque*this.inverseMomentOfInertia*dt;
	this.angularVelocity *= this.ANGULAR_DECAY;
	this.angularPosition += this.angularVelocity*dt;

	this.recomputeModel();
}


Body.prototype.applyImpulse = function(j, n) {
	// j is scalar amount of impulse
	// n is the direction of impulse

	// translational
	this.velocity = this.velocity.plus(n.times(j*this.inverseMass));
	
}

Body.prototype.applyRotationalImpulse = function(j, n, c) {
	// rotational
	var r = c.minus(this.cm);
	var dv = j*this.inverseMomentOfInertia*(r.cross(n));
	this.angularVelocity += dv;
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
	var TOLERANCE = 0;
	var Amin = MAXIMUM_VALUE,
		Amax = MINIMUM_VALUE,
		Bmin = MAXIMUM_VALUE,
		Bmax = MINIMUM_VALUE;
	// n should have been normalized
	n.normalize();

	for (var i = 0; i < A.v.length; ++i) {
		var cur = A.getAbsolutePosition(i).dot(n);
		Amin = Math.min(Amin, cur); 
		Amax = Math.max(Amax, cur);
	}

	if (B instanceof Sphere) {
		var begin = B.cm.plus(n.times(B.radius)).dot(n);
		var end = B.cm.minus(n.times(B.radius)).dot(n);
		Bmin = Math.min(begin, end);
		Bmax = Math.max(begin, end);
	} else {
		for (var i = 0; i < B.v.length; ++i) {
			var cur = B.getAbsolutePosition(i).dot(n);
			Bmin = Math.min(Bmin, cur);
			Bmax = Math.max(Bmax, cur);
		}
	}
	return !(Amax-TOLERANCE <= Bmin || Bmax-TOLERANCE <= Amin);
}

function applyFriction(A, B, n, j, c, vrel) {
	// vrel is relative velocity of B seen from A
	// A is the reference
	var vrel = B.velocity.minus(A.velocity);

	var tangent = vrel.minus(n.times(vrel.dot(n)));
	tangent.normalize();
	//if (Math.abs(tangent.x) < 0.1) tangent.x = 0;

	var rA = c.minus(A.cm);
	var rB = c.minus(B.cm);


	var f = vrel.dot(tangent) / (A.inverseMass + B.inverseMass + A.inverseMomentOfInertia*Math.pow(rA.cross(tangent), 2) + B.inverseMomentOfInertia*Math.pow(rB.cross(tangent), 2));
	
	// combined static friction
	var us = Math.sqrt(A.staticFriction*A.staticFriction+B.staticFriction*B.staticFriction);
	
	var frictionImpulse;
	if (Math.abs(f) < j * us) {
		frictionImpulse = tangent.times(f);
	} else {
		// combined dynamic friction
		var ud = Math.sqrt(A.dynamicFriction*A.dynamicFriction + B.dynamicFriction*B.dynamicFriction);
		frictionImpulse = tangent.times(j * ud);
	}
	A.velocity = A.velocity.plus(frictionImpulse.times(A.inverseMass));
	B.velocity = B.velocity.minus(frictionImpulse.times(B.inverseMass));

	A.angularVelocity -= A.inverseMomentOfInertia*(rA.cross(frictionImpulse));
	B.angularVelocity -= B.inverseMomentOfInertia*(rB.cross(frictionImpulse));
}


function resolveCollision(A, B) {
	// A is the reference
	// if A or B is polygon, swap such that A is polygon
	if (A instanceof Sphere) {
		var C = A;
		A = B;
		B = C;
	}

	// choose a normal from A
	// n must point out of A
	var res = A.axisOfLeastSeparationWith(B);
	var n = res.axis;
	var depth = res.depth;

	var contacts = [];
	if (!(A instanceof Sphere) && !(B instanceof Sphere)) {
		// getting contact of incident B
		var ref = A.getBestContactEdge(n);
		var inc = B.getBestContactEdge(n.flip());
		contacts = getContactPoints(ref, inc);
	} else {
		contacts = B.getContact(A);
		if (!(A instanceof Sphere)) {
			n = res.impactAxis;
		}
	}


	// relative velocity of B as seen from A
	var vrel = B.velocity.minus(A.velocity);
	var proj = vrel.dot(n);

	// for (var i = 0; i < contacts.length; ++i) {	

	// 	var c = contacts[i];
	// 	var rA = contacts[i].minus(A.cm).dot(n);
	// 	var rB = contacts[i].minus(B.cm).dot(n);

	// 	var e = Math.min(A.e, B.e);
	// 	var j = -(1+e)*proj/(A.inverseMass + B.inverseMass + rA*rA*A.inverseMomentOfInertia + rB*rB*B.inverseMomentOfInertia);
	// 	if(proj < 0) {
	// 		A.applyImpulse(-j/contacts.length, n);
	// 		B.applyImpulse(j/contacts.length, n);
	// 		A.applyRotationalImpulse(-j/contacts.length, n, c);
	// 		B.applyRotationalImpulse(j/contacts.length, n, c);
	// 	}
	// 	applyFriction(A, B, n, j/contacts.length);
	// }
	positionalCorrection(A, B, n, depth);

	if (proj > 0) return;

	for (var i = 0; i < contacts.length; ++i) {	
		//var vrel = B.velocity.minus(A.velocity);
		//var proj = vrel.dot(n);
		// if (proj > 1e-9) break;

		var c = contacts[i];
		var rA = contacts[i].minus(A.cm).cross(n);
		var rB = contacts[i].minus(B.cm).cross(n);


		drawContactPoint(c);
		var e = Math.min(A.e, B.e);
		var j = -(1+e)*proj/(A.inverseMass + B.inverseMass + rA*rA*A.inverseMomentOfInertia + rB*rB*B.inverseMomentOfInertia);

		A.applyImpulse(-j, n);
		B.applyImpulse(j, n);
		A.applyRotationalImpulse(-j, n, c);
		B.applyRotationalImpulse(j, n, c);
		applyFriction(A, B, n, j, c, vrel);
	}	

}

function positionalCorrection(A, B, n, depth) {
	// A is the reference body
	// n points out of A
	var CORRECTION_RATE = 1;
	var correction = CORRECTION_RATE * depth/(A.inverseMass + B.inverseMass);
	// if (Math.abs(n.x*correction)<2) n.x = 0;
	A.cm = A.cm.plus(n.flip().times(correction*A.inverseMass));
	B.cm = B.cm.plus(n.times(correction*B.inverseMass));

}

function resolveCollisionExtended(A, B) {
	var resA = A.axisOfLeastSeparationWith(B);
	var resB = B.axisOfLeastSeparationWith(A);
	if (resA.depth < resB.depth) {
		resolveCollision(A, B);
	} else {
	 	resolveCollision(B, A);
	}
}


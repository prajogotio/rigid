function Sphere(cm, radius, 
			  mass, e, momentOfInertia,
			  velocity, acceleration,
			  angularVelocity) {
	this.radius = radius;
	
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

	this.staticFriction = 0.43;
	this.dynamicFriction = 0.33;

	this.ANGULAR_DECAY = 0.90;
}

Sphere.prototype = Object.create(Body.prototype);
Sphere.prototype.constructor = Sphere;

Sphere.prototype.axisOfLeastSeparationWith = function(b) {
	if (b instanceof Sphere) {
		var n = b.cm.minus(this.cm);
		var depth = b.radius + this.radius - n.length();
		n.normalize();
		return {
			axis: n,
			depth: depth
		}
	} else {
		var MAXIMUM_VALUE = 1e301;
		var chosen = -1;
		var contact = null;
		var lowest = MAXIMUM_VALUE;
		for (var i = 0; i < b.v.length; ++i) {
			var j = (i+1)%b.v.length;
			var u = b.getAbsolutePosition(b.v[i]);
			var v = b.getAbsolutePosition(b.v[j]);
			var n = b.getAbsoluteDirection(b.n[i]);
			var res = rayWithSegmentIntersection(this.cm, n, u, v);
			var c;
			if (res.t <= 0) {
				c = u;
			} else if (res.t >= 1) {
				c = v;
			} else {
				c = res.intersect;
			}
			var dist = c.minus(this.cm).length();
			if (dist < lowest) {
				contact = c;
				lowest = dist;
				chosen = i;
			}
		}

		var axis = b.getAbsoluteDirection(b.n[chosen]);
		var depth =  this.radius+(halfPlaneContains(contact, axis, this.cm) ? -1 : 1)*lowest;

		return {
			axis : axis,
			depth : depth,
			index: chosen,
			contact: contact,
			impactAxis : this.cm.minus(contact).normalize(),
		}
	}
}

Sphere.prototype.getContact = function(b) {
	var res = this.axisOfLeastSeparationWith(b);
	var n = res.axis;
	var d = res.depth;
	if (b instanceof Sphere) {
		// R r
		// s t, s+t = d
		// s = (R^2-r^2+d^2)/2d
		var length = this.cm.minus(b.cm).length();
		var s = (this.radius*this.radius-b.radius*b.radius+length*length)/(2*length);
		return [this.cm.plus(n.times(s))];
	} else {
		return [this.axisOfLeastSeparationWith(b).contact];
	}
}

Sphere.prototype.collidesWith = function(b) {
	if (b instanceof Sphere) {
		return this.cm.minus(b.cm).length() <= this.radius + b.radius;
	} else {
		for (var i = 0; i < b.v.length; ++i) {
			var j = (i+1)%b.v.length;
			var u = b.getAbsolutePosition(b.v[i]);
			var v = b.getAbsolutePosition(b.v[j]);
			var n = b.getAbsoluteDirection(b.n[i]);
			var res = rayWithSegmentIntersection(this.cm, n, u, v);
			var c;
			if (res.t <= 0) {
				c = u;
			} else if (res.t >= 1) {
				c = v;
			} else {
				c = res.intersect;
			}
			if (c.minus(this.cm).length() <= this.radius) {
				return true;
			}
		}
		return false;
	}
}
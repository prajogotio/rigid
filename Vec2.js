function Vec2(x, y) {
	this.x = x;
	this.y = y;
};

Vec2.prototype.minus = function(v) {
	return new Vec2(this.x-v.x, this.y-v.y);
}

Vec2.prototype.normalize = function() {
	var length = this.length();
	if (length == 0) return this;
	this.x /= length;
	this.y /= length;
	return this;
}

Vec2.prototype.length = function() {
	return Math.sqrt(this.x*this.x+this.y*this.y);
}

Vec2.prototype.rotate = function(theta) {
	theta = Math.PI/180 * theta;
	var x = this.x;
	var y = this.y;
	var cos = Math.cos(theta);
	var sin = Math.sin(theta);
	this.x = x * cos - y * sin;
	this.y = x * sin + y * cos;
	return this;
}

Vec2.prototype.dot = function(v) {
	return this.x * v.x + this.y * v.y;
}

Vec2.prototype.copy = function() {
	return new Vec2(this.x, this.y);
}

Vec2.prototype.flip = function() {
	return new Vec2(-this.x, -this.y);
}

Vec2.prototype.plus = function(v) {
	return new Vec2(this.x + v.x, this.y + v.y);
}

Vec2.prototype.times = function(c) {
	return new Vec2(this.x * c, this.y * c);
}

Vec2.prototype.cross = function(v) {
	return this.x*v.y - this.y*v.x;
}

function getNormal(u, v) {
	var z = v.minus(u);
	z.rotate(90);
	z.normalize();
	return z;
}

function halfPlaneContains(A, n, B) {
	// is B contained in the half-plane defined by point A and normal n?
	var v = B.minus(A);
	return n.dot(v) >= 0;
}

function rayWithSegmentIntersection(A, n, U, V) {
	// ray at A with direction n
	// segment UV
	var UV = V.minus(U);
	var LHS = n.y * (U.x - A.x) + n.x * (A.y - U.y);
	var RHS = UV.y * n.x - UV.x * n.y;
	if (RHS == 0) {
		if (LHS == 0) {
			return {consistent:true, finite:false};
		} else {
			return {consistent:false};
		}
	} else {
		return {consistent:true, finite:true, t: LHS/RHS, intersect:U.plus(UV.times(LHS/RHS))}
	}
}


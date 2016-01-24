

function getManifolds(A, B) {
	var resA = A.axisOfLeastSeparationWith(B);
	var resB = B.axisOfLeastSeparationWith(A);
	var n = resA.axis;
	var depth = resA.depth;
	if (resA.depth > resB.depth) {
		var C = A;
		A = B;
		B = C;
		n = resB.axis;
		depth = resB.depth;
	} 
	var c = A.getContact(B);
	var manifolds = [];
	for (var i = 0; i < c.length; ++i){
		manifolds.push({'inc': A, 'ref': B, 'contact':c[i], 'axis':n, 'accN':0, 'accT':0});
	}
	positionalCorrection(A, B, n, depth);
	return manifolds;
}

function resolveManifold(m) {
	var A = m.inc;
	var B = m.ref;
	var c = m.contact;
	var n = m.axis;

	var vrel = computeRelativeVelocity(A, B, c);
	var proj = vrel.dot(n);

	var rA = c.minus(A.cm).cross(n);
	var rB = c.minus(B.cm).cross(n);

	var e = Math.min(A.e, B.e);
	var j = -(1+e)*proj/(A.inverseMass + B.inverseMass + rA*rA*A.inverseMomentOfInertia + rB*rB*B.inverseMomentOfInertia);

	var temp = Math.max(m.accN+j, 0);
	var jn = temp - m.accN;

	A.applyImpulse(-jn, n);
	B.applyImpulse(jn, n);
	A.applyRotationalImpulse(-jn, n, c);
	B.applyRotationalImpulse(jn, n, c);
	applyManifoldFriction(m);
	m.accN = temp;
}



function applyManifoldFriction(m) {
	var A = m.inc;
	var B = m.ref;
	var c = m.contact;
	var n = m.axis;
	// vrel is relative velocity of B seen from A
	// A is the reference
	var vrel = computeRelativeVelocity(A, B, c);
	var tangent = vrel.minus(n.times(vrel.dot(n)));
	tangent.normalize();
	//if (Math.abs(tangent.x) < 0.1) tangent.x = 0;
	var rA = c.minus(A.cm);
	var rB = c.minus(B.cm);

	var f = -vrel.dot(tangent) / (A.inverseMass + B.inverseMass + A.inverseMomentOfInertia*Math.pow(rA.cross(tangent), 2) + B.inverseMomentOfInertia*Math.pow(rB.cross(tangent), 2));
	// combined dynamic friction
	var ud = Math.sqrt(A.dynamicFriction*A.dynamicFriction + B.dynamicFriction*B.dynamicFriction);
	
	var sum = Math.min(Math.max(f+m.accT, -ud*m.accN), ud*m.accN);
	var frictionImpulse = sum-m.accT;
	m.accT = sum;

	A.applyImpulse(-frictionImpulse, tangent);
	B.applyImpulse(frictionImpulse, tangent);

	A.applyRotationalImpulse(-frictionImpulse, tangent, c);
	B.applyRotationalImpulse(frictionImpulse, tangent, c);
}

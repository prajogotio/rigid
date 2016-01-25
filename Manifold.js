

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
		manifolds.push({'inc': A, 'ref': B, 'contact':c[i], 'axis':n, depth:depth, 'accN':0, 'accT':0});
	}
	
	return manifolds;
}

function prepareManifold(m) {
	var baumgarte = 0.2;
	var dt = 1/30;
	var slop = 0.01;

	var A = m.inc;
	var B = m.ref;
	var c = m.contact;
	var n = m.axis;

	var rA = c.minus(A.cm).cross(n);
	var rB = c.minus(B.cm).cross(n);
	var e = 0; Math.min(A.e, B.e);
	var massNormal = A.inverseMass + B.inverseMass + rA*rA*A.inverseMomentOfInertia + rB*rB*B.inverseMomentOfInertia;
	var vbias = baumgarte/dt*Math.max(0, m.depth-slop);

	var t = n.copy().rotate(90);
	var rAt = c.minus(A.cm).cross(t);
	var rBt = c.minus(B.cm).cross(t);
	var massTangent = A.inverseMass + B.inverseMass + rAt*rAt*A.inverseMomentOfInertia + rBt*rBt*B.inverseMomentOfInertia;

	m.massNormal = massNormal;
	m.massTangent = massTangent;
	m.vbias = vbias;
	m.e = e;


	A.applyImpulse(-m.accN, n);
	B.applyImpulse(m.accN, n);
	A.applyRotationalImpulse(-m.accN, n, c);
	B.applyRotationalImpulse(m.accN, n, c);

	A.applyImpulse(-m.accT, t);
	B.applyImpulse(m.accT, t);
	A.applyRotationalImpulse(-m.accT, t, c);
	B.applyRotationalImpulse(m.accT, t, c);


}

function resolveManifold(m) {
	var A = m.inc;
	var B = m.ref;
	var c = m.contact;
	var n = m.axis;

	drawContactPoint(c);

	var vrel = computeRelativeVelocity(A, B, c);
	var proj = vrel.dot(n);

	var j = (-(1+m.e)*proj+m.vbias)/m.massNormal;
	var temp = Math.max(m.accN+j, 0);
	var jn = temp - m.accN;
	m.accN = temp;

	A.applyImpulse(-jn, n);
	B.applyImpulse(jn, n);
	A.applyRotationalImpulse(-jn, n, c);
	B.applyRotationalImpulse(jn, n, c);
	applyManifoldFriction(m, jn);

}



function applyManifoldFriction(m, jn) {
	var A = m.inc;
	var B = m.ref;
	var c = m.contact;
	var n = m.axis;
	// vrel is relative velocity of B seen from A
	// A is the reference
	var vrel = computeRelativeVelocity(A, B, c);
	var tangent = n.copy().rotate(90);
	tangent.normalize();
	var rA = c.minus(A.cm);
	var rB = c.minus(B.cm);

	var f = -vrel.dot(tangent) / m.massTangent;
	// combined dynamic friction
	var ud = Math.sqrt(A.dynamicFriction*A.dynamicFriction + B.dynamicFriction*B.dynamicFriction);
	
	var sum = Math.min(Math.max(f+m.accT, -ud*m.accN), ud*m.accN);
	var frictionImpulse = sum - m.accT;
	m.accT = sum;

	A.applyImpulse(-frictionImpulse, tangent);
	B.applyImpulse(frictionImpulse, tangent);

	A.applyRotationalImpulse(-frictionImpulse, tangent, c);
	B.applyRotationalImpulse(frictionImpulse, tangent, c);
}

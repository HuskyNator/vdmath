module vdmath.quaternions;
import vdmath.mat;
import std.math;

alias Quat = Quaternion!float;

struct Quaternion(Type = float) {
	alias QuatType = typeof(this);
	Type r = 1;
	Vec!(3, Type) v = Vec!(3, Type)(0);

	this(Type r, Vec!(3, Type) v) {
		this.r = r;
		this.v = v;
	}

	this(Type r, Type x, Type y, Type z) {
		this.r = r;
		this.v = Vec!(3, Type)([x, y, z]);
	}

	// Conjugate
	QuatType opUnary(string op)() const if (op == "~") {
		return QuatType(r, -v);
	}

	// Inverse
	QuatType opUnary(string op)() const if (op == "-") {
		Type factor = Type(1.0) / (r * r + v.dot(v));
		return QuatType(factor * r, -factor * v);
	}

	QuatType opBinary(string op)(const QuatType r) const if (op == "*") {
		alias l = this;
		return QuatType(l.r * r.r - l.v.dot(r.v), r.v * l.r + l.v * r.r + l.v.cross(r.v));
	}

	// Rotation
	Vec!(3, Type) opBinary(string op)(const Vec!(3, Type) r) const if (op == "^") {
		QuatType qr = QuatType(0, r);
		return (this * qr * ~this).v; // TODO: double check, should be inverse not conjugate, but matters not for non-scaling quaternions. Correct or problematic?
	}

	void opOpAssign(string op, T)(T value) {
		this = mixin(this, op, value.stringof);
	}

	static QuatType rotation(Vec!(3, Type) axis, Type angle) {
		if (angle == 0)
			return QuatType(cast(Type) 1, Vec!(3, Type)(cast(Type) 0));
		return QuatType(cos(angle / 2.0), axis * cast(Type) sin(angle / 2.0));
	}

	auto toMat(uint n = 3)() const if (n == 3 || n == 4) {
		Mat!(n, Type) res = Mat!(n, Type)(1);
		res[0][0] = 1 - 2 * v.y * v.y - 2 * v.z * v.z;
		res[0][1] = 2 * v.x * v.y - 2 * r * v.z;
		res[0][2] = 2 * v.x * v.z + 2 * r * v.y;
		res[1][0] = 2 * v.x * v.y + 2 * r * v.z;
		res[1][1] = 1 - 2 * v.x * v.x - 2 * v.z * v.z;
		res[1][2] = 2 * v.y * v.z - 2 * r * v.x;
		res[2][0] = 2 * v.x * v.z - 2 * r * v.y;
		res[2][1] = 2 * v.y * v.z + 2 * r * v.x;
		res[2][2] = 1 - 2 * v.x * v.x - 2 * v.y * v.y;
		return res;
	}
}

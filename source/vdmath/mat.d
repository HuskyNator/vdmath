module vdmath.mat;
import vdmath.misc;

import std.conv : to;
import std.exception : enforce;
import std.math : abs, cos, sin, sqrt;
import std.traits : isCallable, isFloatingPoint, ReturnType;

alias Vec(size_t size = 3, Type = float) = Mat!(size, 1, Type);
alias Mat(size_t dimension = 3, Type = float) = Mat!(dimension, dimension, Type);

struct Mat(size_t row_count, size_t column_count, Type = float)
		if (row_count > 0 && column_count > 0) {
	enum ulong size = row_count * column_count; // #elements, not #bytes!
	enum bool isVec = (column_count == 1);
	enum bool isMat = !isVec;
	enum bool isSquare = (column_count == row_count);
	enum bool isArithmetic = __traits(isArithmetic, Type);

	alias MatType = typeof(this);

	union {
		static if (isArithmetic)
			Type[size] vec = 0; // Standard arithmetic values are 0
		else
			Type[size] vec;
		Type[column_count][row_count] mat;
		static if (isVec) {
			struct {
				static if (size >= 1)
					Type x;
				static if (size >= 2)
					Type y;
				static if (size >= 3)
					Type z;
				static if (size >= 4)
					Type w;
			}
		}
	}

	static if (isVec)
		alias this = vec;
	else
		alias this = mat;

	this(Type val) {
		static if (isVec)
			vec[] = val;
		else {
			foreach (i; 0 .. row_count)
				foreach (j; 0 .. column_count)
					static if (isArithmetic)
						mat[i][j] = (i == j) ? val : 0;
					else
						mat[i][j] = (i == j) ? val : Type.init;
		}
	}

	this(Type[size] vals...) {
		this.vec = vals;
	}

	this(Type[column_count][row_count] vals...) {
		this.mat = vals;
	}

	unittest {
		Vec!3 v1 = Vec!3(2.0f);
		assert(v1.vec == [2.0f, 2.0f, 2.0f]);
		Mat!3 m1 = Mat!3(2.0f);
		assert(m1.mat == [
				[2.0f, 0.0f, 0.0f], [0.0f, 2.0f, 0.0f], [0.0f, 0.0f, 2.0f]
			]);

		Vec!3 v2 = Vec!3(1.0f, 2.0f, 3.0f);
		assert(v2.vec == [1.0f, 2.0f, 3.0f]);
		Mat!2 m2 = Mat!2(1.0f, 2.0f, 3.0f, 4.0f);
		assert(m2.mat == [[1.0f, 2.0f], [3.0f, 4.0f]]);

		float[2][2] vals4M = [[1.0f, 2.0f], [3.0f, 4.0f]];
		float[4] vals4V = [1.0f, 2.0f, 3.0f, 4.0f];
		Mat!2 m3 = Mat!2(vals4M);
		assert(m3.vec == vals4V);
		assert(m3.mat == vals4M);
	}

	static if (isVec)
		this(size_t L, size_t R)(Type[L] left, Type[R] right...) if (L + R == size) {
			this.vec = left ~ right;
		}

	unittest {
		Vec!2 a = Vec!2(1, 2);
		Vec!3 b = Vec!3(3, 4, 5);
		Vec!4 c = Vec!4(a, 3, 4);
		assert(c == [1, 2, 3, 4]);
		Vec!5 d = Vec!5(a, b);
		assert(d == [1, 2, 3, 4, 5]);

	}

	void setCol(size_t k, Vec!(row_count, Type) col) {
		assert(k < row_count);
		foreach (i; 0 .. row_count) {
			this.mat[k][i] = col[i];
		}
	}

	Vec!(row_count, Type) col(size_t index) {
		Vec!(row_count, Type) column;
		foreach (r; 0 .. row_count)
			column.vec[r] = this.mat[r][index];
		return column;
	}

	Vec!(column_count, Type) row(size_t index) {
		return Vec!(column_count, Type)(this.mat[index]);
	}

	static if (isSquare && isArithmetic) {
		alias ArithType = ResultType!(Type, "+", Type); // Integer Promotion

		// TODO: add determinant function?

		Mat!(row_count, column_count, ArithType) inverse() {
			Mat!(row_count, column_count, ArithType) inverse;
			ArithType determinant;
			static if (row_count == 2) {
				inverse[0][0] = mat[1][1];
				inverse[0][1] = -mat[0][1];
				inverse[1][0] = -mat[1][0];
				inverse[1][1] = mat[0][0];
				determinant = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
			} else static if (row_count == 3) {
				// Determine adjugate (transposed cofactor) matrix (2x2 determinants times even index sign)
				inverse[0][0] = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
				inverse[0][1] = -(mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]);
				inverse[0][2] = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];

				inverse[1][0] = -(mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]);
				inverse[1][1] = mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0];
				inverse[1][2] = -(mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]);

				inverse[2][0] = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];
				inverse[2][1] = -(mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]);
				inverse[2][2] = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];

				determinant = mat[0][0] * inverse[0][0] + mat[0][1] * inverse[1][0] + mat[0][2] * inverse[2][0];
			} else static if (row_count == 4) {
				// Determine 2x2 determinants for bottom 3 rows
				// Mij_kl references the top left index ij & bottom right index kl
				ArithType M10_21 = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];
				ArithType M10_31 = mat[1][0] * mat[3][1] - mat[1][1] * mat[3][0];
				ArithType M20_31 = mat[2][0] * mat[3][1] - mat[2][1] * mat[3][0];

				ArithType M10_22 = mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0];
				ArithType M10_32 = mat[1][0] * mat[3][2] - mat[1][2] * mat[3][0];
				ArithType M20_32 = mat[2][0] * mat[3][2] - mat[2][2] * mat[3][0];

				ArithType M10_23 = mat[1][0] * mat[2][3] - mat[1][3] * mat[2][0];
				ArithType M10_33 = mat[1][0] * mat[3][3] - mat[1][3] * mat[3][0];
				ArithType M20_33 = mat[2][0] * mat[3][3] - mat[2][3] * mat[3][0];

				ArithType M11_22 = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
				ArithType M11_32 = mat[1][1] * mat[3][2] - mat[1][2] * mat[3][1];
				ArithType M21_32 = mat[2][1] * mat[3][2] - mat[2][2] * mat[3][1];

				ArithType M11_23 = mat[1][1] * mat[2][3] - mat[1][3] * mat[2][1];
				ArithType M11_33 = mat[1][1] * mat[3][3] - mat[1][3] * mat[3][1];
				ArithType M21_33 = mat[2][1] * mat[3][3] - mat[2][3] * mat[3][1];

				ArithType M12_23 = mat[1][2] * mat[2][3] - mat[1][3] * mat[2][2];
				ArithType M12_33 = mat[1][2] * mat[3][3] - mat[1][3] * mat[3][2];
				ArithType M22_33 = mat[2][2] * mat[3][3] - mat[2][3] * mat[3][2];

				// Determine adjugate (transposed cofactor) matrix (minor times even index sign)
				// Using a laplace expansion to determine the minor
				inverse[0][0] = mat[1][1] * M22_33 - mat[1][2] * M21_33 + mat[1][3] * M21_32;
				inverse[0][1] = -(mat[0][1] * M22_33 - mat[0][2] * M21_33 + mat[0][3] * M21_32);
				inverse[0][2] = mat[0][1] * M12_33 - mat[0][2] * M11_33 + mat[0][3] * M11_32;
				inverse[0][3] = -(mat[0][1] * M12_23 - mat[0][2] * M11_23 + mat[0][3] * M11_22);

				inverse[1][0] = -(mat[1][0] * M22_33 - mat[1][2] * M20_33 + mat[1][3] * M20_32);
				inverse[1][1] = mat[0][0] * M22_33 - mat[0][2] * M20_33 + mat[0][3] * M20_32;
				inverse[1][2] = -(mat[0][0] * M12_33 - mat[0][2] * M10_33 + mat[0][3] * M10_32);
				inverse[1][3] = mat[0][0] * M12_23 - mat[0][2] * M10_23 + mat[0][3] * M10_22;

				inverse[2][0] = mat[1][0] * M21_33 - mat[1][1] * M20_33 + mat[1][3] * M20_31;
				inverse[2][1] = -(mat[0][0] * M21_33 - mat[0][1] * M20_33 + mat[0][3] * M20_31);
				inverse[2][2] = mat[0][0] * M11_33 - mat[0][1] * M10_33 + mat[0][3] * M10_31;
				inverse[2][3] = -(mat[0][0] * M11_23 - mat[0][1] * M10_23 + mat[0][3] * M10_21);

				inverse[3][0] = -(mat[1][0] * M21_32 - mat[1][1] * M20_32 + mat[1][2] * M20_31);
				inverse[3][1] = mat[0][0] * M21_32 - mat[0][1] * M20_32 + mat[0][2] * M20_31;
				inverse[3][2] = -(mat[0][0] * M11_32 - mat[0][1] * M10_32 + mat[0][2] * M10_31);
				inverse[3][3] = mat[0][0] * M11_22 - mat[0][1] * M10_22 + mat[0][2] * M10_21;

				// determine determinant with cofactors
				determinant = mat[0][0] * inverse[0][0] + mat[0][1] * inverse[1][0] + mat[0][2]
					* inverse[2][0] + mat[0][3] * inverse[3][0];
			}

			enforce(determinant != 0, "Matrix not invertable: determinant = 0");
			inverse = inverse / determinant;
			return inverse;
		}

		unittest {
			Mat!2 identity2 = Mat!2(1);
			Mat!3 identity3 = Mat!3(1);
			Mat!4 identity4 = Mat!4(1);

			Mat!2 m1 = Mat!2([3, 5, 7, 11]);
			Mat!2 i1 = m1.inverse().mult(m1);
			assert(i1.almostEquals(identity2));

			Mat!3 m2 = Mat!3([1, 0, 3, 4, 5, 6, 7, 8, 9]); // Determinant -12
			Mat!3 i2 = m2.inverse().mult(m2);
			assert(i2.almostEquals(identity3));

			Mat!4 m3 = Mat!4([
				1, -2, 3, 4, 5, 6, 7, -8, 9, 10, 11, 12, 13, 14, 15, 16
			]); // Determinant 512
			Mat!4 i3 = m3.inverse().mult(m3);
			assert(i3.almostEquals(identity4));
		}
	}

	Mat!(column_count, row_count, Type) transposed() const {
		Mat!(column_count, row_count, Type) result;
		foreach (i; 0 .. row_count)
			foreach (j; 0 .. column_count)
				result.mat[j][i] = this.mat[i][j];
		return result;
	}

	unittest {
		Mat!(3, 2) m1 = Mat!(3, 2)([[1.0f, 2.0f], [3.0f, 4.0f], [5.0f, 6.0f]]);
		Mat!(2, 3) correct = Mat!(2, 3)([[1.0f, 3.0f, 5.0f], [2.0f, 4.0f, 6.0f]]);
		assert(m1.transposed() == correct);
	}

	// programming language problem: I cant define non-auto result type, because its ResultType!(Type, "+", other.Type),
	// but we dont have other (yet), nor would we have access to other.Type, apart from defining the type of other to have a generic T
	// would be nice to be able to be able to simply query its Type "alias/template parameter"
	// Would also be nice to be able to define result type after having access to other. (does Rust do this? given return type is after)
	// Would it be possible to predefine a type and define its actual value after?
	// Should be possible to do with explicit templates no that i think about it. . . . ? ! ? ! ? . . .

	template dot(R : T[size], T) if (isVec) {
		alias ResultT = ResultType!(Type, "*", T);
		ResultT dot(const R right) const {
			ResultT result = 0;
			static foreach (i; 0 .. size)
				result += this.vec[i] * right[i];
			return result;
		}
	}

	unittest {
		Vec!4 v1 = Vec!4(1, 2, 3, 4);
		Vec!4 v2 = Vec!4(2, 3, 4, 5);
		float[4] arr2 = [2, 3, 4, 5];
		float correct = 40;
		assert(v1.dot(v2) == correct);
		assert(v1.dot(arr2) == correct);
	}

	template cross(R : T[size], T) if (isVec && row_count == 3) {
		alias ArithType = ResultType!(Type, "*", T);
		alias ResultT = Vec!(3, ArithType);

		ResultT cross(const R right) const {
			ResultT result;
			result.vec[0] = this.vec[1] * right[2] - right[1] * this.vec[2];
			result.vec[1] = this.vec[2] * right[0] - right[2] * this.vec[0];
			result.vec[2] = this.vec[0] * right[1] - right[0] * this.vec[1];
			return result;
		}
	}

	unittest {
		Vec!3 v1 = Vec!3(1, 2, 3);
		Vec!3 v2 = Vec!3(2, 3, 4);
		Vec!3 correct = Vec!3(-1, 2, -1);
		assert(v1.cross(v2) == correct);
	}

	static if (isVec) {
		static if (is(Type == real))
			alias LengthT = real;
		else static if (is(Type == double))
			alias LengthT = double;
		else
			alias LengthT = float;

		T length2(T = LengthT)() const {
			T len = 0;
			foreach (i; 0 .. size)
				len += this.vec[i] ^^ 2;
			return len;
		}

		unittest {
			Vec!3 v1 = Vec!3(1, 2, 3);
			assert(v1.length2() == 14);
		}

		T length(T = LengthT)() const {
			return sqrt(length2!T());
		}

		unittest {
			Vec!3 v1 = Vec!3(1, 2, 2);
			assert(v1.length() == 3);
			Vec!3 v2 = Vec!3(1, 2, 3);
			assert(abs(v2.length() - 3.74165738677) < 1e-5);
		}

		Vec!(size, T) normalize(T = LengthT)() const {
			Vec!(size, T) n;
			n.vec[] = this.vec[];
			n.vec[] = n.vec[] * cast(T)(1.0 / this.length());
			return n;
		}

		unittest {
			Vec!3 v1 = Vec!3(1, 2, 2);
			Vec!3 correct = Vec!3(1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0);
			v1.normalize().assertAlmostEquals(correct);
		}
	}

	template mult(T, size_t K) {
		alias ArithType = ResultType!(Type, "*", T);
		alias ResultT = Mat!(row_count, K, ArithType);

		ResultT mult(const T[K][column_count] right) const {
			ResultT result = 0;
			foreach (i; 0 .. row_count)
				foreach (j; 0 .. K)
					foreach (k; 0 .. column_count)
						result.mat[i][j] += this.mat[i][k] * right[k][j];
			return result;
		}
	}

	auto mult(T)(const T[column_count] right) const
	if (is(ResultType!(Type, "*", T))) {
		return mult!(T, 1)(cast(T[1][column_count]) right);
	}

	unittest {
		Mat!(2, 3) m1 = Mat!(2, 3)([[1, 2, 3], [2, 3, 4]]);
		Mat!(3, 2) m2 = Mat!(3, 2)([1, 2], [3, 4], [5, 6]);
		Mat!(2, 2) correct1 = Mat!(2, 2)([22, 28], [31, 40]);
		assert(m1.mult(m2) == correct1);

		Vec!3 v1 = Vec!3(1, 2, 3);
		Vec!2 correct2 = Vec!2(14, 20);
		assert(m1.mult(v1) == correct2);
	}

	MatType opUnary(string op)() const if (op == "-") {
		return this * (cast(Type)-1);
	}

	unittest {
		Vec!3 v1 = Vec!3(1, 2, 3);
		Vec!3 correct = Vec!3(-1, -2, -3);
		assert(-v1 == correct);
	}

	auto opBinary(string op, R)(const R right) const if (op == "^") {
		return this.mult(right);
	}

	unittest {
		Mat!(2, 3) m1 = Mat!(2, 3)([[1, 2, 3], [2, 3, 4]]);
		Vec!3 v1 = Vec!3(1, 2, 3);
		Vec!2 correct1 = Vec!2(14, 20);
		assert((m1 ^ v1) == correct1);
	}

	auto opOpAssign(string op, T)(T value) {
		this = mixin(this, op, "value");
		return this;
	}

	unittest {
		Mat!2 m1 = Mat!2(1, 2, 3, 4);
		m1 += 1;
		Mat!2 correct = Mat!2(2, 3, 4, 5);
		assert(m1 == correct);
	}

	T sum(T = float)() const {
		T result = 0;
		static foreach (i; 0 .. size)
			result += this.vec[i];
		return result;
	}

	unittest {
		Mat!2 a;
		foreach (i; 0 .. a.size)
			a.vec[i] = 2 * i;
		assert(a.sum() == 12);
	}

	static if (__traits(compiles, Type.max))
		Type min() const {
			Type minVal = Type.max;
			foreach (i; 0 .. size)
				if (vec[i] < minVal)
					minVal = vec[i];
			return minVal;
		}

	unittest {
		Vec!3 a = Vec!3(-2, 1, 2);
		assert(a.min() == -2);
	}

	static if (__traits(isFloating, Type) || __traits(compiles, Type.min))
		Type max() const {
			static if (__traits(isFloating, Type))
				Type maxVal = -Type.max; // floating points don't have a 'min' value and are signed magnitude.
			else
				Type maxVal = Type.min;
			foreach (i; 0 .. size)
				if (vec[i] > maxVal)
					maxVal = vec[i];
			return maxVal;
		}

	unittest {
		Vec!3 a = Vec!3(-2, 1, 2);
		assert(a.max() == 2);
	}

	// OpApply provided by alias this.

	unittest {
		int sum = 0;
		void count(float add) {
			sum += cast(int) add;
		}

		Vec!3 v1 = Vec!3(1, 2, 3);
		foreach (v; v1)
			count(v);
		int correct = 6;
		assert(sum == correct);
	}

	// TODO: re-evaluate "dual-context" deprecation
	auto map(alias fun)() const
	if (isCallable!fun && __traits(compiles, fun(Type.init))) {
		Mat!(row_count, column_count, ReturnType!fun) result;
		foreach (i; 0 .. size)
			result.vec[i] = fun(this.vec[i]);
		return result;
	}

	unittest {
		bool even(float number) {
			return (cast(int) number) % 2 == 0;
		}

		Mat!2 m1 = Mat!2(1, 2, 3, 4);
		auto answer = m1.map!even();
		Mat!(2, bool) correct = Mat!(2, bool)(false, true, false, true);

		assert(is(typeof(answer) == Mat!(2, bool)));
		assert(answer == correct);
	}

	pure nothrow @nogc auto opBinary(string op, R)(const R right) const
	if (is(ResultType!(Type, op, R))) {
		alias ResT = ResultType!(Type, op, R);
		Mat!(row_count, column_count, ResT) result;
		mixin("result.vec[] = this.vec[] " ~ op ~ " right;");
		return result;
	}

	unittest {
		Mat!2 m1 = Mat!2(1, 2, 3, 4);
		Mat!2 result = (m1 + 1) * 2;
		Mat!2 correct = Mat!2(4, 6, 8, 10);
		assert(result == correct);
	}

	static if (isVec) {
		pure nothrow @nogc auto opBinary(string op, RT)(const RT[size] right) const
				if (!isList!(RT) && is(ResultType!(Type, op, RT))) { // WARNING: why cant size be used directly?
			alias ResT = ResultType!(Type, op, RT);
			Mat!(row_count, column_count, ResT) result;
			foreach (i; 0 .. size)
				mixin("result.vec[i] = this.vec[i] " ~ op ~ " right[i];");
			return result;
		}

		unittest {
			Vec!3 v1 = Vec!3(1, 1, 1);
			float[3] v2 = [0.0f, 1.0f, 2.0f];
			Vec!3 correct = Vec!3(1, 2, 3);
			assert(v1 + v2 == correct); // Vec + []

			Vec!3 v3 = Vec!3(v2);
			assert(v1 + v3 == correct); // Vec + Vec
		}
	} else {
		pure nothrow @nogc auto opBinary(string op, R:
			RT[column_count][row_count], RT)(const R right) const
		if (is(ResultType!(Type, op, RT)) && !isList!(RT)) {
			alias ResT = ResultType!(Type, op, RT);
			Mat!(row_count, column_count, ResT) result;
			foreach (i; 0 .. row_count)
				foreach (j; 0 .. column_count)
					mixin("result.mat[i][j] = this.mat[i][j] " ~ op ~ " right[i][j];");
			return result;
		}

		unittest {
			Mat!2 m1 = Mat!2(1, 1, 1, 1);
			float[2][2] m2 = [[0.0f, 1.0f], [2.0f, 3.0f]];
			Mat!2 correct = Mat!2(1, 2, 3, 4);
			assert(m1 + m2 == correct); // Mat + [][]

			Mat!2 m3 = Mat!2(m2);
			assert(m1 + m3 == correct); // Mat + Mat
		}
	}

	M opCast(M : Mat!(row_count, column_count, T), T)() const {
		M result;
		foreach (i; 0 .. size)
			result.vec[i] = this.vec[i].to!T;
		return result;
	}

	unittest {
		Vec!2 v1 = Vec!2(1, 2);
		Vec!(2, int) v2 = cast(Vec!(2, int)) v1;
		string[2] v3 = cast(Vec!(2, string)) v2; // implict alias this
		string[2] correct = ["1", "2"];
		assert(v3 == correct);
	}

	// TODO: revise toString functionality

	static if (isMat)
		string toString(bool nice = false)() const {
			char[] cs;
			cs.reserve(6 * size);
			cs ~= '{';
			static foreach (i; 0 .. row_count) {
				cs ~= '[';
				static foreach (j; 0 .. column_count) {
					cs ~= this.mat[i][j].to!string;
					static if (j != column_count - 1)
						cs ~= ", ";
				}
				static if (i != row_count - 1)
					cs ~= nice ? "],\n " : "], ";
				else
					cs ~= ']';
			}
			cs ~= '}';
			return cast(string) cs;
		}

	static if (isVec)
		string toString() const {
			char[] cs = ['['];
			foreach (v; vec)
				cs ~= v.to!string ~ ", ";
			cs = cs[0 .. $ - 2] ~ ']';
			return cast(string) cs;
		}

	static if (isVec) {
		bool opEquals(RT)(const RT[size] other) const @safe pure nothrow {
			foreach (size_t i; 0 .. size)
				if (this.vec[i] != other[i])
					return false;
			return true;
		}

		unittest {
			Vec!3 v1 = Vec!3(1, 2, 3);
			int[3] v2 = [1, 2, 3];
			int[3] v3 = [1, 2, 2];
			assert(v1 == v1);
			assert(v1 == v2);
			assert(v1 != v3);
		}
	} else {
		// bool opEquals(RT)(const Mat!(row_count, column_count, RT) other) const @safe pure nothrow {
		// 	return this.opEquals(other.vec);
		// }

		// TODO: re-evaluate bug
		// May have template deduction bug: https://forum.dlang.org/post/rjnywrpcsipgkronwrrc@forum.dlang.org
		bool opEquals(RT)(const RT[column_count][row_count] other) const @safe pure nothrow {
			return Vec!(size, Type)(this.vec).opEquals(cast(RT[size]) other);
		}

		unittest {
			Mat!2 m1 = Mat!2(1, 2, 3, 4);
			int[2][2] m2 = [[1, 2], [3, 4]];
			int[2][2] m3 = [[1, 2], [3, 3]];
			assert(m1 == m1);
			assert(m1 == m2);
			assert(m1 != m3);
		}
	}

	static if (isVec) {
		bool almostEquals(RT)(const RT[size] right, double delta = 1e-5) const
		if (__traits(compiles, abs(Type.init - RT.init))) {
			alias SumType = ResultType!(Type, "-", RT);
			SumType diff = 0;
			foreach (i; 0 .. size)
				diff += abs(this.vec[i] - right[i]);
			return diff < delta;
		}

		unittest {
			Vec!4 v1 = Vec!4(1.25f, 2.25f, 2.75f, 3.75f);
			int[4] v2 = [1, 2, 3, 4];
			assert(v1.almostEquals(v2, 1.1f));
			assert(!v1.almostEquals(v2, 1.0f));
		}
	} else {
		bool almostEquals(RT)(const RT[column_count][row_count] right, double delta = 1e-5) const
				if (__traits(compiles, abs(Type.init - RT.init))) {
			return Vec!(size, Type)(this.vec).almostEquals(cast(RT[size]) right, delta);
		}

		unittest {
			Mat!2 m1 = Mat!2(1.25f, 2.25f, 2.75f, 3.75f);
			int[2][2] m2 = [[1, 2], [3, 4]];
			assert(m1.almostEquals(m2, 1.1f));
			assert(!m1.almostEquals(m2, 1.0f));
		}
	}

	static if (isVec) {
		void assertAlmostEquals(RT)(const RT[size] right, double delta = 1e-5) const
				if (__traits(compiles, this.almostEquals(right, delta))) {
			alias SumType = ResultType!(Type, "-", RT);
			SumType diff = 0;
			foreach (i; 0 .. size)
				diff += abs(this.vec[i] - right[i]);
			bool eq = diff < delta;
			assert(eq, "Expected " ~ this.toString ~ " ==(delta=" ~ delta.to!string ~ ") " ~ to!string(
					right) ~ " but found difference: " ~ diff.to!string ~ ")");
		}

		unittest {
			import std.exception : assertThrown, assertNotThrown;
			import core.exception : AssertError;

			Vec!4 v1 = Vec!4(1.25f, 2.25f, 2.75f, 3.75f);
			int[4] v2 = [1, 2, 3, 4];

			assertNotThrown!AssertError(v1.almostEquals(v2, 1.1f));
			assertThrown!AssertError(v1.assertAlmostEquals(v2, 1.0f));
		}
	} else {
		void assertAlmostEquals(RT)(const RT[column_count][row_count] right, double delta = 1e-5) const
				if (__traits(compiles, this.almostEquals(right, delta))) {
			Vec!(size, Type)(vec).assertAlmostEquals(cast(RT[size]) right, delta);
		}

		unittest {
			import std.exception : assertThrown, assertNotThrown;
			import core.exception : AssertError;

			Mat!2 m1 = Mat!2(1.25f, 2.25f, 2.75f, 3.75f);
			int[2][2] m2 = [[1, 2], [3, 4]];

			assertNotThrown!AssertError(m1.assertAlmostEquals(m2, 1.1f));
			assertThrown!AssertError(m1.assertAlmostEquals(m2, 1.0f));
		}
	}

	// TODO: re-evaluate hashes
	// Hashes required for associative lists
	static if (is(Type == byte) || is(Type == ubyte) || is(Type == short) || is(Type == ushort)
		|| is(Type == int) || is(Type == uint) || is(Type == long) || is(Type == ulong)) {
		size_t toHash() const @safe pure nothrow {
			size_t hash = 1;
			foreach (Type s; this.vec)
				hash = 31 * hash + s;
			return hash;
		}
	}

	static if (is(Type == bool)) {
		size_t toHash() const @safe pure nothrow {
			size_t hash = 1;
			foreach (Type s; this.vec)
				hash = 31 * hash + s ? 5 : 3;
			return hash;
		}
	}

	static if (is(Type == float)) {
		private static int _castFloatInt(const float f) @trusted {
			return *cast(int*)&f;
		}

		size_t toHash() const @safe pure nothrow {
			size_t hash = 1;
			foreach (Type s; this.vec)
				hash = 31 * hash + _castFloatInt(s); // Reinterpret as int
			return hash;
		}
	}

	static if (is(Type == double)) {
		private static long _castDoubleLong(const double d) @trusted {
			return *cast(long*)&d;
		}

		size_t toHash() const @safe pure nothrow {
			size_t hash = 1;
			foreach (Type s; this.vec)
				hash = 31 * hash + _castDoubleLong(s); // Reinterpret as long
			return hash;
		}
	}

	static if (isSquare && (row_count == 3 || row_count == 4) && isFloatingPoint!Type) {
		static MatType rotationMx(Type angle) {
			MatType rotationM = MatType(1);
			Type cos = cos(angle);
			Type sin = sin(angle);
			rotationM[1][1] = cos;
			rotationM[1][2] = -sin;
			rotationM[2][1] = sin;
			rotationM[2][2] = cos;
			return rotationM;
		}

		unittest {
			import std.math : PI, PI_2;

			Mat!4 rotation = Mat!(4).rotationMx(0);
			Mat!4 correct = Mat!4(1);
			assert(rotation == correct);

			rotation = Mat!(4).rotationMx(PI_2);
			correct = Mat!4(0);
			correct[0][0] = 1;
			correct[1][2] = -1;
			correct[2][1] = 1;
			correct[3][3] = 1;
			assert(rotation.almostEquals(correct));

			rotation = Mat!(4).rotationMx(PI);
			correct = Mat!4(1);
			correct[1][1] = -1;
			correct[2][2] = -1;
			assert(rotation.almostEquals(correct));
		}

		static MatType rotationMy(Type angle) {
			MatType rotationM = MatType(1);
			Type cos = cos(angle);
			Type sin = sin(angle);
			rotationM[0][0] = cos;
			rotationM[0][2] = sin;
			rotationM[2][0] = -sin;
			rotationM[2][2] = cos;
			return rotationM;
		}

		unittest {
			import std.math : PI, PI_2;

			Mat!4 rotation = Mat!(4).rotationMy(0);
			Mat!4 correct = Mat!4(1);
			assert(rotation == correct);

			rotation = Mat!(4).rotationMy(PI_2);
			correct = Mat!4(0);
			correct[0][2] = 1;
			correct[1][1] = 1;
			correct[2][0] = -1;
			correct[3][3] = 1;
			assert(rotation.almostEquals(correct));

			rotation = Mat!(4).rotationMy(PI);
			correct = Mat!4(1);
			correct[0][0] = -1;
			correct[2][2] = -1;
			assert(rotation.almostEquals(correct));
		}

		static MatType rotationMz(Type angle) {
			MatType rotationM = MatType(1);
			Type cos = cos(angle);
			Type sin = sin(angle);
			rotationM[0][0] = cos;
			rotationM[0][1] = -sin;
			rotationM[1][0] = sin;
			rotationM[1][1] = cos;
			return rotationM;
		}

		unittest {
			import std.math : PI, PI_2;

			Mat!4 rotation = Mat!(4).rotationMz(0);
			Mat!4 correct = Mat!4(1);
			assert(rotation == correct);

			rotation = Mat!(4).rotationMz(PI_2);
			correct[0][0] = 0;
			correct[0][1] = -1;
			correct[1][0] = 1;
			correct[1][1] = 0;
			assert(rotation.almostEquals(correct));

			rotation = Mat!(4).rotationMz(PI);
			correct = Mat!4(1);
			correct[0][0] = -1;
			correct[1][1] = -1;
			assert(rotation.almostEquals(correct));

		}
	}

	// TODO: re-evaluate
	// Get rotation required to rotate a vector from the y axis towards the direction.
	// Rotation about the x axis is applied before rotation about z
	// Thus there is no rotation around the y axis.
	// _
	static Vec!3 getRotation(Vec!3 direction) {
		import std.math : acos, atan, PI_2, signbit;

		if (direction.x == 0 && direction.y == 0) {
			if (direction.z == 0)
				return direction; // [0,0,0] -> [0,0,0]
			return Vec!3([PI_2, 0, 0]); // [0,0,z] -> [PI/2,0,0]
		}
		double R = Vec!2([direction.x, direction.y]).length!double();
		// [x,y,z] -> [atan(z/sqrt(x²+y²)),0,-teken(x)*acos(y/sqrt(x²+y²))]
		return Vec!3([
			atan(direction.z / R), 0, -signbit(direction.x) * acos(direction.y / R)
		]);
	}

	// TODO: add test

	// TODO: re-evaluate
	// Opposite of getRotation
	// ROtation around y axis assumed 0
	static Vec!3 getDirection(Vec!3 rotation) {
		import std.math : sin, cos;

		return Vec!3([-sin(rotation.z), cos(rotation.z), sin(rotation.x)]);
	}

	// TODO: add test
}

// Instantiate Templates for Unittests
unittest {
	import std.meta : AliasSeq;

	static foreach (i; 1 .. 10) {
		static foreach (type; AliasSeq!(bool, int, float, double, string)) {
			mixin("Vec!(" ~ i.stringof ~ ',' ~ type.stringof ~ ") v" ~ i.stringof ~ type.stringof ~ ';');
			mixin("Mat!(" ~ i.stringof ~ ',' ~ type.stringof ~ ") m" ~ i.stringof ~ type.stringof ~ ';');
		}
	}

	// TODO: rewrite previous tests to use Type (instead of eg. only float/int)
}

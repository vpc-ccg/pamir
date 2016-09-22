/// 786

#pragma once

#include <fmt/format.h>

#define EN(msg,...)\
	fmt::print(stderr, msg, ##__VA_ARGS__)

#define E(msg,...)\
	fmt::print(stderr, msg"\n", ##__VA_ARGS__)

#define assertit(a,f,...) if (!(a)) {\
	E("Assertion failed: {} at {}:{} as: " f, #a, __FILE__, __LINE__, ##__VA_ARGS__); \
}

#define asserteq(a,b,f,...) if ((a) != (b)) {\
	E("Assertion failed: {} ({}) != {} ({}) at {}:{} as: " f, #a, (a), #b, (b), __FILE__, __LINE__, ##__VA_ARGS__); \
}

inline char getDNAValue (char ch) {
	#define _ 0
	static char c[128] = {
		_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,// 0 15
		_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,// 16 31
		_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,// 32 47
		_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,// 48 63
		_,0,_,1,_,_,_,2,_,_,_,_,_,_,_,_,// 64 79
		_,_,_,_,3,_,_,_,_,_,_,_,_,_,_,_,// 80 95
		_,0,_,1,_,_,_,2,_,_,_,_,_,_,_,_,// 96 111
		_,_,_,_,3,_,_,_,_,_,_,_,_,_,_,_,// 112 127
	};	
	#undef _
	return c[ch];
}

# Big-int
Signed big integer (BigInteger java class analogue) implementation on C++

# Details of implementation

The number is stored as a vector of unsigned ints which are basically the bits of the number in two's complement.

The implementation includes: 
- default constructor
- copy constructor
- constructors from other numeric types
- explicit constructor from std::string
- assignment operator
- comparison operators
- arithmetic operations: addition, subtraction, division, unary minus and plus
- increments and decrements
- bit operations: and, or, xor, not (the behaviour is similar to int type)
- bit shift
- to_string method

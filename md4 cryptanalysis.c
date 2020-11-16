/**
 * Implantation of Hans Dobbertin's cryptanalyis of MD4 (FSE 1996)
 * For further information related to MD4, see rfc1320
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#define F(x, y, z) (((x) & (y)) | ((~x) & (z)))
#define G(x, y, z) (((x) & (y)) | ((x) & (z)) | ((y) & (z)))
#define H(x, y, z) ((x) ^ (y) ^ (z))
#define ROTATE_LEFT(x, n) (((x) << (n)) | ((x & 0xFFFFFFFF) >> (32-(n))))
#define J(a,b,c,d,e,f)  (F(a, b, 0) - F(c, d, -1) - ROTATE_LEFT(e, 13) + ROTATE_LEFT(f, 13))

#define FF(a, b, c, d, x, s) {a = ROTATE_LEFT((a + x + F(b,c,d)),s);}
#define GG(a, b, c, d, x, s) {a = ROTATE_LEFT((a + x + G(b,c,d) + (unsigned long int)0x5a827999),s);}
#define HH(a, b, c, d, x, s) {a = ROTATE_LEFT((a + x + H(b,c,d) + (unsigned long int)0x6ed9eba1),s);}

/* used during continuous approximation step */
unsigned int distance(unsigned int j){
	int i,d=0;
	for(i=0;i<32;i++)
		d += ((j & (1<<i))>>i);

	return d;
}

void md4_dobbertin(){
  srand(getpid());

  int B11  =  0;
  int A12  =  -1;
  int A12_ =  0;
  
  int C14, B15, A16, D17, C18, B19;
  int R14, R15, R16, R17, R18, R19;
  int B15_, C14_, D13, D13_;
  int R15_, R14_, R13, R13_;

  int A5, B5, C5, D5, t;

  int A8, D9, C10;
  int X0, X1, X2,  X3,  X4,  X5,  X6,  X7,
	  X8, X9, X10, X11, X12, X13, X14, X15;

  unsigned int j;

  int timer;
  int a1,b1,c1,d1;
  int a2,b2,c2,d2;
	
step1:

  timer = 0;

  do{

	C14 = rand();
    B15 = rand();
    A16 = rand();
    D17 = rand();
    C18 = rand();
    B19 = rand();

    B15_ = B15 - G(C18 + 32, D17, A16) + G(C18, D17, A16) + ROTATE_LEFT(B19 - 33554432, 19) - ROTATE_LEFT(B19, 19) - 1;
    C14_ = C14 - G(D17, A16, B15_) + G(D17, A16, B15) + ROTATE_LEFT(C18 + 32, 23) - ROTATE_LEFT(C18, 23);
    D13 = ROTATE_LEFT(C14, 21) - ROTATE_LEFT(C14_, 21);
    D13_ = D13 - G(A16, B15_, C14_) + G(A16, B15, C14);

  }
  while( (G(B15, C14, D13) - G(B15_, C14_, D13_)) != 1);

	R13 = D13;	R13_ = D13_;
  	R14 = C14;	R14_ = C14_;
  	R15 = B15;	R15_ = B15_;
  	R16 = A16;
  	R17 = D17;
  	R18 = C18;
  	R19 = B19;

	j = J(R14_,R13_,R14,R13,R15_,R15);

	goto step2;

step2:	/* continuous approximation step */

  while(distance(j) > 1){

	if(timer > 500)	/* the best upper bound for timer has been found empirically */
		goto step1;

	R14 = C14 ^ (1 << (rand() & 31));
	R15 = B15 ^ (1 << (rand() & 31));
	R16 = A16 ^ (1 << (rand() & 31));
	R17 = D17 ^ (1 << (rand() & 31));
	R18 = C18 ^ (1 << (rand() & 31));
	R19 = B19 ^ (1 << (rand() & 31));

    R15_ = R15 - G(R18 + 32, R17, R16) + G(R18, R17, R16) + ROTATE_LEFT(R19 - 33554432, 19) - ROTATE_LEFT(R19, 19) - 1;
    R14_ = R14 - G(R17, R16, R15_) + G(R17, R16, R15) + ROTATE_LEFT(R18 + 32, 23) - ROTATE_LEFT(R18, 23);
    R13 =  ROTATE_LEFT(R14, 21) - ROTATE_LEFT(R14_, 21);
    R13_ = R13 - G(R16, R15_, R14_) + G(R16, R15, R14);

	if( 
	(distance(J(R14_,R13_,R14,R13,R15_,R15)) < distance(j))	/* we didn't manage to use original version of continuous approximation, so we modified it */
	&&
	((G(R15, R14, R13) - G(R15_, R14_, R13_)) == 1)
	){

	C14 = R14;	C14_ = R14_;
	B15 = R15;	B15_ = R15_;
	A16 = R16;
	D17 = R17;
	C18 = R18;
	B19 = R19;
	D13 = R13;	D13_ = R13_;

	j = J(C14_,D13_,C14,D13,B15_,B15);

	}

	timer++;
  }  

  goto step3;

step3:

	C14 -= j;
	B15_ = B15 - G(C18 + 32, D17, A16) + G(C18, D17, A16) + ROTATE_LEFT(B19 - 33554432, 19) - ROTATE_LEFT(B19, 19) - 1;
    C14_ = C14 - G(D17, A16, B15_) + G(D17, A16, B15) + ROTATE_LEFT(C18 + 32, 23) - ROTATE_LEFT(C18, 23);
    D13 = ROTATE_LEFT(C14, 21) - ROTATE_LEFT(C14_, 21);
    D13_ = D13 - G(A16, B15_, C14_) + G(A16, B15, C14);

	j = J(C14_,D13_,C14,D13,B15_,B15);

	if(j || ((G(B15, C14, D13) - G(B15_, C14_, D13_)) != 1))
		goto step1;

	goto step4;

step4:

	if(G(B19, C18, D17) != G(B19 - 33554432, C18 + 32, D17))
	goto step1;

	printf(" inner almost collision found \n");

	goto step5;

step5:

	X13 = 0x74697069;
	C10 = ROTATE_LEFT(D13_, 25) - ROTATE_LEFT(D13, 25);
	X14 = ROTATE_LEFT(C14, 21) - C10 - F(D13, -1, 0);
	X15 = ROTATE_LEFT(B15, 13) - B11 - F(C14, D13, -1);
	X0  = ROTATE_LEFT(A16, 29) + 1   - G(B15, C14, D13) - 0x5a827999;
	X4  = ROTATE_LEFT(D17, 27) - D13 - G(A16, B15, C14) - 0x5a827999;
	X8  = ROTATE_LEFT(C18, 23) - C14 - G(D17, A16, B15) - 0x5a827999;
	X12 = ROTATE_LEFT(B19, 19) - B15 - G(C18, D17, A16) - 0x5a827999;
	D9  = ROTATE_LEFT(D13, 25) - F(-1, 0, C10) - X13;
	A8  = ROTATE_LEFT( -1, 29) - F(B11, C10, D9) - X12;

	goto step6;

step6:	/* differential attack */

	X1  = rand();
	X2  = rand();
	X3  = rand();
	X5  = rand();

	A5 = 0x67452301;
	B5 = 0xefcdab89; 
	C5 = 0x98badcfe;
	D5 = 0x10325476;

	FF(A5, B5, C5, D5, X0,  3);
	FF(D5, A5, B5, C5, X1,  7);
	FF(C5, D5, A5, B5, X2, 11);
	FF(B5, C5, D5, A5, X3, 19);
	FF(A5, B5, C5, D5, X4,  3);
	FF(D5, A5, B5, C5, X5,  7);

	t   = ROTATE_LEFT(A8,  29) - A5 - X8;
	X6  = ROTATE_LEFT(t,   21) - C5 - F(D5, A5, B5);
	X7  = -1 - B5 - F(t, D5, A5);
	X9  = ROTATE_LEFT(D9,  25) - D5 - F(A8, -1, t);
	X10 = ROTATE_LEFT(C10, 21) - t - F(D9, A8, -1);
	X11 = 1 - F(C10, D9, A8);

	a1 = A16; b1 = B19; c1 = C18; d1 = D17;
	a2 = A16; b2 = B19 - 33554432; c2 = C18 +32; d2 = D17;

    GG (a1, b1, c1, d1, X1,   3); /* 20 */
    GG (d1, a1, b1, c1, X5,   5); /* 21 */
    GG (c1, d1, a1, b1, X9,   9); /* 22 */
	GG (b1, c1, d1, a1, X13, 13); /* 23 */
	GG (a1, b1, c1, d1, X2,   3); /* 24 */
	GG (d1, a1, b1, c1, X6,   5); /* 25 */
	GG (c1, d1, a1, b1, X10,  9); /* 26 */
	GG (b1, c1, d1, a1, X14, 13); /* 27 */
	GG (a1, b1, c1, d1, X3,   3); /* 28 */
	GG (d1, a1, b1, c1, X7,   5); /* 29 */
	GG (c1, d1, a1, b1, X11,  9); /* 30 */
	GG (b1, c1, d1, a1, X15, 13); /* 31 */
	HH (a1, b1, c1, d1, X0,   3); /* 32 */
	HH (d1, a1, b1, c1, X8,   9); /* 33 */
	HH (c1, d1, a1, b1, X4,  11); /* 34 */
	HH (b1, c1, d1, a1, X12, 15); /* 35 */

	GG (a2, b2, c2, d2, X1,   3); /* 20 */
    GG (d2, a2, b2, c2, X5,   5); /* 21 */
    GG (c2, d2, a2, b2, X9,   9); /* 22 */
	GG (b2, c2, d2, a2, X13, 13); /* 23 */
	GG (a2, b2, c2, d2, X2,   3); /* 24 */
	GG (d2, a2, b2, c2, X6,   5); /* 25 */
	GG (c2, d2, a2, b2, X10,  9); /* 26 */
	GG (b2, c2, d2, a2, X14 ,13); /* 27 */
	GG (a2, b2, c2, d2, X3,   3); /* 28 */
	GG (d2, a2, b2, c2, X7,   5); /* 29 */
	GG (c2, d2, a2, b2, X11,  9); /* 30 */
	GG (b2, c2, d2, a2, X15, 13); /* 31 */
	HH (a2, b2, c2, d2, X0,   3); /* 32 */
	HH (d2, a2, b2, c2, X8,   9); /* 33 */
	HH (c2, d2, a2, b2, X4,  11); /* 34 */
	HH (b2, c2, d2, a2, X12 +1, 15); /* 35 */
						
	if(
		(a1 == a2)&&
		(b1 == b2)&&
		(c1 == c2)&&
		(d1 == d2)
		){

	  printf(" collision found \n");
	  printf(" Message 1 = %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x \n \n", X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15);
	  printf(" Message 2 = %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x \n \n", X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12 + 1, X13, X14, X15);

	  goto step6;
	}
	goto step6;
	
}

int main(int argc, char **argv){
	md4_dobbertin();
	exit(0);
}

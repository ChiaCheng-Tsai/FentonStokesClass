
#ifndef FentonHeaders_H
#define FentonHeaders_H

#define Skip(stream)	fgets(dummy,400,stream)
#define Readtext(stream,x)	fgets(x,400,stream);    x[strlen(x)-1] = '\0'
#define Read(stream,x,y) {fscanf(stream,"%"#y, &x);Skip(stream);}
#define Fenton_pi			3.14159265358979324
#define iff(x,y)	if(strcmp(x,#y)==0)
#define LO			"\t%7.4f" /* "\t%10.6f" */

#endif


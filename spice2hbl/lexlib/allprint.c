#include <stdio.h>
int printable(int);

void allprint(c)
  char c; {
	extern FILE *yyout;
	switch(c){
		case '\n':
			fprintf(yyout,"\\n");
			break;
		case '\t':
			fprintf(yyout,"\\t");
			break;
		case '\b':
			fprintf(yyout,"\\b");
			break;
		case ' ':
			fprintf(yyout," ");
			break;
		default:
			if(!printable(c))
				fprintf(yyout,"\\%-3o",c);
			else
				putc(c,yyout);
			break;
		}
	return;
	}
void sprint(s)
  char *s; {
	while(*s)
		allprint(*s++);
	return;
	}
int printable(c)
  int c;
	{
	return((040 < c && c < 0177) || (c&0300)==0300 );
	}

# include <stdio.h>
allprint(c)
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
sprint(s)
  char *s; {
	while(*s)
		allprint(*s++);
	return;
	}
printable(c)
  int c;
	{
	return((040 < c && c < 0177) || (c&0300)==0300 );
	}

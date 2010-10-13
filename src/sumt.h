int      DoCompareTree (void);
int      DoCompareTreeParm (char *parmName, char *tkn);
int      DoSumt (void);
int      DoSumtParm (char *parmName, char *tkn);
int      DoSumtTree (void);
int      DoSumtTreeParm (char *parmName, char *tkn);
void     ResetTranslateTable (void);
int		 ShowConTree (FILE *fp, PolyTree *t, int screenWidth, int showSupport);
void     ShowParts (FILE *fp, SafeLong *p, int nTaxaToShow);

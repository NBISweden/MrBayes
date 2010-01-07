

extern int      AddToString (char *s1, char *s2, int *x);
extern unsignedSafeLong Expecting (int y);
extern int      CheckString (char *s1, char *s2, int *x);
extern int      CheckStringValidity (char *s);
int             DoExecute (void);
extern void     FinishTree (TreeNode *p, int *i, int isThisTreeRooted);
extern int      FreeMatrix (void);
extern void     GetUserDownPass (TreeNode *p, TreeNode **x, int *y);
void     		GetToken (int *tokenType);
int      		FindValidCommand (char *tk, int *numMatches);
extern int      IsArgValid (char *s, char *validArg);
int             IsIn (char ch, char *s);
int             IsSame (char *s1, char *s2);
int             IsWhite (char c);
extern int      ParseCommand (char *s);
extern int      RootUserTree (TreeNode *p);
extern void     SetUpParms (void);
void     		ShowNodes (TreeNode *p, int indent, int isThisTreeRooted);
int      		ShowTree (TreeNode *r, int isThisTreeRooted, int nTips);
void            State_CODON (char *state, int code, int division);
void            State_DOUBLET (char *state, int code);
int				StateCode_AA (int n);
int				StateCode_NUC4 (int n);
int				StateCode_Std (int n);
char     		WhichNuc (int x);

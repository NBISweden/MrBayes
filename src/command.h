#ifndef __COMMAND_H__
#define __COMMAND_H__

int         AddString (char ***list, int len, char *token);
BitsLong    Expecting (int y);
int         CheckString (char **list, int len, char *token, int *matchIndex);
int         CheckStringValidity (char *s);
int         DoExecute (void);
int         FreeMatrix (void);
int         GetToken (char *token, int *tokenType, char **sourceH);
int         FindValidCommand (char *tk, int *numMatches);
int         IsArgValid (char *s, char *validArg);
int         IsIn (char ch, char *s);
int         IsSame (char *s1, char *s2);
int         IsWhite (char c);
int         ParseCommand (char *s);
void        ResetCharacterFlags (void);
void        ResetTaxaFlags (void);
int         RootUserTree (TreeNode *p);
void        SetUpParms (void);
void        ShowNodes (TreeNode *p, int indent, int isThisTreeRooted);
int         ShowTree (Tree *t);
void        State_CODON (char *state, int code, int division);
void        State_DOUBLET (char *state, int code);
int         StateCode_AA (int n);
int         StateCode_NUC4 (int n);
int         StateCode_Std (int n);
char        WhichAA (int x);
char        WhichNuc (int x);
char        WhichRes (int x);
char        WhichStand (int x);

#endif  /* __COMMAND_H__ */

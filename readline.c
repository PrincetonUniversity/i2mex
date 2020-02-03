/*  dummy libreadline.a routines -- when gnu readline & command  */
/*  editing are not provided for command line interactive codes  */
 
char *readline(char *prompt)
 
{
  return (char *) 1;    /* this signals the dummy routine */
}
 
void add_history()
{
  return;
}
 
void history_list()
{
  return;
}

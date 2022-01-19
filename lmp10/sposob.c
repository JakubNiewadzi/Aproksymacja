int spr(int x, double **tab){
  int pom;
  for (int i = 0; i < x; i++ )
  {
    pom=0;
    for (int j = 0; j < x; j++ )
    {
      if(tab[i][j]!=0)
      {
        pom=1;
      }
    }
    if(pom==0)
    {
      
      return 1;
    }
  }

  for (int i = 0; i < x; i++ )
  {
    pom=0;
    for (int j = 0; j < x; j++ )
    {
      if(tab[j][i]!=0)
      {
        pom=1;
      }
    }
    if(pom==0)
    {

      return 1;
    }
  }
  return 0;
}

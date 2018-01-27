loadData[fname_]:=Module[{etab, et},
			If[FileExistsQ[fname] && (FileByteCount[fname]>0),
			   etab=Select[Import[fname,{"Table"}],#[[3]]=!=0.&];
			   If[etab[[-1]][[1]]=="Not",etab=etab[[1;;-2]]];
			   et=etab;
			   et=Transpose[etab];
			   et[[3]]=et[[3]]-et[[3]][[1]];
			   et=Transpose[et];			   
			   et]
		       ];

getFreeEnergyFromData[et_,t_]:=Module[{ Z, Ev, cv},
			     Z=Plus@@Map[Exp[#[[3]]-#[[1]]/t]&,et];
			     Ev=(Plus@@Map[#[[1]]*Exp[#[[3]]-#[[1]]/t]&,et])/Z;
			     cv=((Plus@@Map[(#[[1]]^2)*Exp[#[[3]]-#[[1]]/t]&,et])/Z-Ev^2)/t^2;
			     {Log[Z],Ev,cv}]

FileExistsQ::fstr: 
   File specification fname is not a string of one or more characters.

FileByteCount::fstr: 
   File specification fname is not a string of one or more characters.

   
getFreeEnergy[fname_,t_]:=Module[{etab, Z, Ev, cv, et},
				 If[FileExistsQ[fname] && (FileByteCount[fname]>0),
				    etab=Select[Import[fname,{"Table"}],#[[3]]=!=0.&];
				    If[etab[[-1]][[1]]=="Not",etab=etab[[1;;-2]]];
				    et=etab[[2;;]];
(*				    et=Transpose[etab];
				    et[[3]]=et[[3]]-et[[3]][[1]];
				    et=Transpose[et];*)
				    Z=Plus@@Map[Exp[#[[3]]-#[[1]]/t]&,et];
				    Ev=(Plus@@Map[#[[1]]*Exp[#[[3]]-#[[1]]/t]&,et])/Z;
				    cv=((Plus@@Map[(#[[1]]^2)*Exp[#[[3]]-#[[1]]/t]&,et])/Z-Ev^2)/t^2;
				    {Log[Z],Ev,cv},
				    0]
				];

Z[M_,N_,K_,t_]:=Product[(1-Exp[-(i+j+k-1)/t])/(1-Exp[-(i+j+k-2)/t]),{i,M},{j,N},{k,K}]

Z[M_,N_,K_]:=Product[(1-Exp[-(i+j+k-1)/t])/(1-Exp[-(i+j+k-2)/t]),{i,M},{j,N},{k,K}]


Z2[M_,N_,K_]:=Product[(i+j+k-1)/(i+j+k-2),{i,M},{j,N},{k,K}]

freeEn[M_,N1_,K_]:=N[1/(M*N1)*Log[Z2[M,N1,K]]]

Out[128]= 0.748933

freeEn[10,10,2]

Export["wang-math-plot-30.png", 
Plot[freeEn[10,x,10],{x,10,10000},PlotRange->{{0,10000},{0,1}}]
      ]

                    
Out[135]= wang-math-plot-31.png

Export["wang-math-plot-31.png", 
Plot[freeEn[20,x,20],{x,10,10000},PlotRange->{All,All}]
      ]

2+2

[Calculating...]

Out[137]= 4

Out[136]= 4


                    
Out[134]= wang-math-plot-30.png

                    
Out[133]= wang-math-plot-30.png

                    
Out[132]= wang-math-plot-30.png

                    
Out[131]= wang-math-plot-30.png

                    
Out[130]= wang-math-plot-30.png

Out[129]= 0.225023

          Log[20]
Out[126]= -------
             4

Export["wang-math-plot-3.png", 
Plot[t*Log[Z[5,5,5,t]],{t,0.05,5}]
      ]

                  
Out[50]= wang-math-plot-3.png
 
Export["wang-math-plot-4.png", 
Plot[q*D[q*Log[Z[5,5,5,{q,2}]],q]/.{q->t},{t,0.05,15}]
]


expr=Simplify[q*D[q*Log[Z[5,5,5,{q,2}]],q]/.{q->t}]


q*D[q*Log[Z[5,5,5]],{q,2}]

Out[65]= 0

expr=q*D[q*Log[Z[5,5,5,q]],{q,2}]/.{q->t};

expr=t*D[t*Log[Z[30,30,30,t]],{t,2}];

expr3=t*D[t*Log[Z[10,10,10,t]],{t,2}];


Export["wang-math-plot-12.png", 
Plot[expr3,{t,0.05,15}]
      ]

2+2

Out[80]= 4



                  
Out[79]= wang-math-plot-12.png
 
                  
Out[77]= wang-math-plot-11.png

nterrupt

Out[74]= $Aborted


Export["wang-math-plot-11.png", 
Plot[expr,{t,0.05,15}]
      ]


wang-math-plot-11.png
 

Export["wang-math-plot-8.png", 
Plot[expr,{t,0.05,15}]
      ]

Export["wang-math-plot-10.png", 
Plot[expr,{t,0.05,15}]
      ]

                  
Out[73]= wang-math-plot-10.png
    
expr2=-q^2*D[Log[Z[5,5,5,q]],{q,1}]/.{q->t};
 

Export["wang-math-plot-9.png", 
Plot[expr2,{t,0.05,5}]
]

                  
Out[71]= wang-math-plot-9.png
 
                  
Out[70]= wang-math-plot-9.png
 
                  
Out[69]= wang-math-plot-9.png
 
Out[67]= wang-math-plot-8.png

 
expr

Out[64]= 0
                  
Out[63]= wang-math-plot-8.png
 

expr

                         1/t        2/t        3/t        4/t        5/t
Out[57]= {t ((125 + 124 E    + 119 E    + 234 E    + 462 E    + 554 E    + 
 
                6/t         7/t         8/t         9/t         10/t
>          751 E    + 1070 E    + 1341 E    + 1590 E    + 2047 E     + 
 
                 11/t         12/t         13/t         14/t         15/t
>          2489 E     + 2945 E     + 3302 E     + 3740 E     + 4223 E     + 
 
                 16/t         17/t         18/t         19/t         20/t
>          4654 E     + 4913 E     + 5205 E     + 5484 E     + 5558 E     + 
 
                 21/t         22/t         23/t         24/t         25/t
>          5570 E     + 5585 E     + 5444 E     + 5181 E     + 4915 E     + 
 
                 26/t         27/t         28/t         29/t         30/t
>          4555 E     + 4192 E     + 3766 E     + 3295 E     + 2837 E     + 
 
                 31/t         32/t         33/t         34/t         35/t
>          2471 E     + 2027 E     + 1635 E     + 1323 E     + 1055 E     + 
 
                36/t        37/t        38/t        39/t        40/t
>          761 E     + 578 E     + 410 E     + 284 E     + 180 E     + 
 
                41/t       42/t       43/t       44/t      45/t    46/t
>          124 E     + 71 E     + 38 E     + 16 E     + 6 E     + E    ) / 
 
                1/t        4/t        2/t    4/t
>        ((1 + E   ) (1 + E   ) (1 - E    + E   ) 
 
                 1/t    2/t    3/t    4/t        3/t    6/t
>          (1 - E    + E    - E    + E   ) (1 + E    + E   ) 
 
                 1/t    2/t    3/t    4/t    5/t    6/t
>          (1 - E    + E    - E    + E    - E    + E   ) 
 
                 1/t    2/t    3/t    4/t    5/t    6/t    7/t    8/t
>          (1 + E    + E    + E    + E    + E    + E    + E    + E    + 
 
              9/t    10/t        1/t    2/t    3/t    4/t    5/t    6/t
>            E    + E    ) (1 + E    + E    + E    + E    + E    + E    + 
 
              7/t    8/t    9/t    10/t    11/t    12/t
>            E    + E    + E    + E     + E     + E    ) t) + 
 
                   -14/t         -13/t 2       -12/t 3        -11/t 4
>       Log[((1 - E     ) (-1 + E     )  (1 - E     )  (-1 + E     )  
 
                  -10/t 5       -9/t 3       -8/t
>           (1 - E     )  (1 - E    )  (1 - E    )) / 
 
                 -7/t        -6/t 3       -5/t 5        -4/t 4       -3/t 3
>         ((1 - E    ) (1 - E    )  (1 - E    )  (-1 + E    )  (1 - E    )  
 
                   -2/t 2       -(1/t)
>           (-1 + E    )  (1 - E      ))]), 
 
                  -7         -(13/2) 2       -6 3        -(11/2) 4       -5 5
>    t Log[((1 - E  ) (-1 + E       )  (1 - E  )  (-1 + E       )  (1 - E  )  
 
                -(9/2) 3       -4   2
>         (1 - E      )  (1 - E  ) E ) / 
 
               -(7/2)        -3 3       -(5/2) 5        -2 4       -(3/2) 3
>       ((1 - E      ) (1 - E  )  (1 - E      )  (-1 + E  )  (1 - E      )  
 
                  1             2
>         (1 - -------) (-1 + E) )]}
               Sqrt[E]

?D

D[Sin[x],{x,2}]

Out[60]= -Sin[x]

Out[59]= 6 x

D[f, x] gives the partial derivative \[PartialD] f/\[PartialD] x
                                                            n                n
    . D[f, {x, n}] gives the multiple derivative \[PartialD]  f/\[PartialD] x

      . D[f, x, y, ...] differentiates f

         successively with respect to x, y, ...
         .D[f, {{x , x , ...}}] for a scalar f
                  1   2
            gives the vector derivative 
            (\[PartialD] f/\[PartialD] x , \[PartialD] f/\[PartialD] x , ...)
                                        1                             2
            . D[f, {array}] gives a tensor derivative.


Export["wang-math-plot-7.png", 
Plot[expr,{t,0.05,15}]
]

                  
Out[56]= wang-math-plot-7.png
 

Out[53]= wang-math-plot-4.png

                  
Out[52]= wang-math-plot-4.png
 
?D

D[f, x] gives the partial derivative \[PartialD] f/\[PartialD] x
                                                            n                n
    . D[f, {x, n}] gives the multiple derivative \[PartialD]  f/\[PartialD] x

      . D[f, x, y, ...] differentiates f

         successively with respect to x, y, ...
         .D[f, {{x , x , ...}}] for a scalar f
                  1   2
            gives the vector derivative 
            (\[PartialD] f/\[PartialD] x , \[PartialD] f/\[PartialD] x , ...)
                                        1                             2
            . D[f, {array}] gives a tensor derivative.

                  
Out[44]= wang-math-plot-4.png
 
                  
Out[39]= wang-math-plot-4.png
 
N[D[t*Log[Z[5,5,5,t]],t]/.{t->1}]

Out[38]= 3.20731
                  
General::ivar: 0.0501011 is not a valid variable.

General::ivar: 0.151122 is not a valid variable.

General::ivar: 0.252142 is not a valid variable.

General::stop: Further output of General::ivar
     will be suppressed during this calculation.

Out[37]= wang-math-plot-4.png

                  
General::ivar: 0.0501011 is not a valid variable.

General::ivar: 0.151122 is not a valid variable.

General::ivar: 0.252142 is not a valid variable.

General::stop: Further output of General::ivar
     will be suppressed during this calculation.

Out[36]= wang-math-plot-4.png


D[t*Log[Z[5,5,5]],t]

D[t*Log[Z[5,5,5,t]],t]

                  
General::ivar: 0.0501011 is not a valid variable.

General::ivar: 0.151122 is not a valid variable.

General::ivar: 0.252142 is not a valid variable.

General::stop: Further output of General::ivar
     will be suppressed during this calculation.

Out[32]= wang-math-plot-4.png
                  
Out[28]= wang-math-plot-3.png

 

dt=loadData["5-5-5-20-3.txt"];


dt2=loadData["10-10-10-20.txt"];

                  
Out[83]= wang-math-plot-13.png
 
Export["wang-math-plot-13.png", 
ListPlot[Table[{t,t*getFreeEnergyFromData[dt2,t][[1]]},{t,0.05,5,0.05}]]
      ]



Export["wang-math-plot-15.png", 
ListPlot[Table[{t,getFreeEnergyFromData[dt2,t][[3]]},{t,0.05,15,0.05}]]
]


dt3=loadData["ising-data-10-10-08-20.txt"];


Export["wang-math-plot-16.png", 
ListPlot[Table[{t,getFreeEnergyFromData[dt3,t][[3]]},{t,0.05,5,0.05}]]
]


dt4=loadData["ising-data-5-5-08-20.txt"];

                  
Out[93]= wang-math-plot-17.png


Export["wang-math-plot-18.png", 
ListPlot[{Table[{t,getFreeEnergyFromData[dt4,t][[3]]},{t,0.05,5,0.05}],Table[{t,getFreeEnergyFromData[dt3,t][[3]]},{t,0.05,5,0.05}]},PlotRange->{All,All}]
      ]

dt6=loadData["noboundary-10-10-10-20.txt"];

expr2=t*D[t*Log[Z[10,10,10,t]],{t,2}];

Export["theory-vs-experimets.png", 
Show[ListPlot[{Table[{t,getFreeEnergyFromData[dt6,t][[3]]},{t,0.05,15,0.05}]},PlotRange->{All,All}],Plot[expr2,{t,0.05,15}]]
      ]

Export["theory-vs-experiment.eps", 
Show[ListPlot[{Table[{t,getFreeEnergyFromData[dt6,t][[3]]},{t,0.05,15,0.05}]},PlotRange->{All,All}],Plot[expr2,{t,0.05,15}]]
      ]

                  
Out[12]= theory-vs-experiment.eps

                  
Out[11]= theory-vs-experiment.pdf


Out[10]= theory-vs-experimets.png


dt5=loadData["noboundary-5-5-5-20.txt"];

                    
Out[115]= wang-math-plot-19.png


expr1=t*D[t*Log[Z[5,5,5,t]],{t,2}];

                    
Out[118]= wang-math-plot-20.png


Export["wang-math-plot-12.png", 
Plot[expr3,{t,0.05,15}]
      ]


Out[122]= 2
                  
Export["wang-math-plot-22.png", 
Show[ListPlot[{Table[{t,getFreeEnergyFromData[dt5,t][[3]]},{t,0.05,5,0.05}]},PlotRange->{All,All}],Plot[expr1,{t,0.05,15}]]
      ]

Out[121]= 2

                    
Out[120]= wang-math-plot-22.png

                    
Out[119]= wang-math-plot-21.png

                    
Out[116]= wang-math-plot-19.png


Ordering[list] gives the positions in list
     at which each successive element of Sort[list]
      appears. Ordering[list, n]

       gives the positions in list
        at which the first n elements of Sort[list]
          appear. Ordering[list, -n]

           gives the positions of the last n
            elements of Sort[list]. Ordering[list, n, p] uses Sort[list, p]. 


Table[{t,getFreeEnergyFromData[dt4,t][[3]]},{t,0.05,5,0.05}]

Ordering[Table[{t,getFreeEnergyFromData[dt3,t][[3]]},{t,0.05,5,0.05}], -1, #[[2]]&]

Out[102]= {5., 9.71903}

ttt=Transpose[Table[{t,getFreeEnergyFromData[dt3,t][[3]]},{t,0.05,5,0.05}]];

Out[112]= {47}

ttt=Transpose[Table[{t,getFreeEnergyFromData[dt4,t][[3]]},{t,0.05,5,0.05}]];

Max[ttt[[2]]]

Ordering[ttt[[2]],-1]

Out[113]= 2.35

ttt[[1]][[47]]

Out[110]= 2.4

Out[109]= 23.1379

Out[108]= {48}

Out[107]= 23.1379

Out[106]= 5.

Table[{t,getFreeEnergyFromData[dt3,t][[3]]},{t,0.05,5,0.05}][[100]]

Max[x , x , ...] yields the numerically largest of the x
     1   2                                              i
    . Max[{x , x , ...}, {y , ...}, ...]
            1   2          1
      yields the largest element of any of the lists. 

Out[103]= {5., 9.71903}

Table[{t,getFreeEnergyFromData[dt4,t][[3]]},{t,0.05,5,0.05}][[100]]

Out[101]= {100}

Out[100]= {5., 2.66294}

?Sort

Out[99]= {100}

Sort[list] sorts the elements of list
     into canonical order. Sort[list, p] sorts using the ordering function p. 

                  
Out[95]= wang-math-plot-18.png

                  
Out[94]= wang-math-plot-18.png


                  
Out[91]= wang-math-plot-16.png

                  
Out[89]= wang-math-plot-15.png

                  
Out[88]= wang-math-plot-15.png

                  
Out[87]= wang-math-plot-15.png

                  
Out[86]= wang-math-plot-14.png

                  
Out[85]= wang-math-plot-14.png
 

                  
Out[84]= wang-math-plot-14.png
 

dt=dt[[2;;]];

getFreeEnergyFromData[dt,10]


Export["wang-math-plot-5.png", 
ListPlot[Table[{t,t*getFreeEnergyFromData[dt,t][[1]]},{t,0.05,5,0.05}]]
]

                  
Out[49]= wang-math-plot-5.png
 
                  
Out[48]= wang-math-plot-5.png
 
                  
Out[45]= wang-math-plot-5.png
 

Export["wang-math-plot-6.png", 
ListPlot[Table[{t,getFreeEnergyFromData[dt,t][[3]]},{t,0.05,5,0.05}]]
]

                  
Out[43]= wang-math-plot-6.png
                  
Out[41]= wang-math-plot-5.png
 
 
Export["wang-math-plot-0.png", 
ListPlot[Table[{t,t*getFreeEnergyFromData[dt,t][[1]]},{t,0.05,5,0.05}]]
      ]

                  
Out[31]= wang-math-plot-0.png
 
                  
Out[30]= wang-math-plot-0.png
 
                  
Out[25]= wang-math-plot-0.png
 

Export["wang-math-plot-2.png", 
ListPlot[Table[{t,getFreeEnergyFromData[dt,t][[3]]},{t,0.05,5,0.05}]]
]


                  
Out[12]= wang-math-plot-2.png
 
                  
Out[11]= wang-math-plot-1.png
                  
Out[10]= wang-math-plot-0.png
 
 

Out[9]= {1197.923001266730, 1.127200936388, 0.0128718026674}

Out[8]= {1197.003451627711, 1.004545793173, 0.007224166193}

Out[7]= {1098.000000000000000, 1.000000000000000, 0.}


Out[6]= {Indeterminate, Indeterminate, Indeterminate}

        
Out[20]= wang-math-plot-4.png

         
Out[16]= wang-math-plot-3.png

         
Out[14]= wang-math-plot-3.png

getFreeEnergy["test.dat",10]


Export["wang-math-plot-0.png", 
ListPlot[Table[{t,getFreeEnergy["test.dat",t][[1]]},{t,0.05,5,0.05}]]
      ]

         
Out[10]= wang-math-plot-1.png
 

Out[8]= -Graphics-

                                                
Set::write: Tag Times in wang-math-plot-0.png -Graphics- is Protected.

Out[9]= -wang-math-plot-0.png
    
 
                               -13
Out[7]= {-30., 300., 2.91038 10   }

Out[6]= {-3000.000000000000, 300.000000000, 0.}

Out[5]= {-300., 300., 0.}



Export["wang-math-plot-1.png", 
ListPlot[Table[{t,getFreeEnergy["test.dat",t][[3]]},{t,0.05,5,0.05}]]]

         
Out[11]= wang-math-plot-2.png
 

Export["wang-math-plot-2.png", 
ListPlot[Table[{t,getFreeEnergy["test.dat",t][[2]]},{t,0.05,5,0.05}]]]

         
Out[12]= wang-math-plot-3.png
 

Export["wang-math-plot-3.png", 
ListPlot[Table[{t,getFreeEnergy["test2.dat",t][[1]]},{t,0.05,5,0.05}]]]

         
Out[18]= wang-math-plot-4.png
 
         
Out[17]= wang-math-plot-3.png


Export["wang-math-plot-4.png", 
ListPlot[Table[{t,getFreeEnergy["test2.dat",t][[3]]},{t,0.05,5,0.05}]]]

dt=loadData["test3.dat"];

                                            -875                  -875
Out[30]= {9377.00000000000, 1.84799633554 10    , 3.94472940532 10    }getFreeEnergyFromData[dt,1]

         
Out[31]= wang-math-plot-5.png
 

Export["wang-math-plot-5.png", 
ListPlot[Table[{t,getFreeEnergyFromData[dt,t][[3]]},{t,0.05,5,0.05}]]]


Export["wang-math-plot-6.png", 
ListPlot[Table[{t,getFreeEnergyFromData[dt,t][[1]]},{t,0.05,5,0.05}]]
      ]

getFreeEnergyFromData[dt,1]/(20*20)

getFreeEnergyFromData[dt,200]/(20*20)

dt2=loadData["test.dat"];

                                              -5
Out[39]= {84650.88358608751, 0.75000000, 0. 10  }

getFreeEnergyFromData[dt2,100]/(20*20)

dt3=loadData["test2.dat"];

                                         -29928280              -29928281
Out[42]= {172296.9548360875, 8.5477152 10         , 6.8381722 10         }

getFreeEnergyFromData[dt3,1000]/(20*20)

[Calculating...]

                                          -29928276              -29928280
Out[44]= {172296.9548360875, 1.14490727 10         , 9.1592582 10         }

                                         -29928624              -29928621
Out[43]= {172296.9548360875, 9.3458838 10         , 7.4767071 10         }

                                              -9
Out[40]= {84651.62608608751, 0.75000000, 0. 10  }

                                             -874                  -877
Out[37]= {23.44250000000000, 1.19050400398 10    , 2.82492529691 10    }

                                             -877                  -877
Out[36]= {23.44250000000000, 1.90728156183 10    , 1.88858447188 10    }

                                             -878                 -878
Out[35]= {23.44250000000000, 4.61999083885 10    , 9.8618235133 10    }



                                             -878                 -878
Out[34]= {23.44250000000000, 4.61999083885 10    , 9.8618235133 10    }

                                            -875                  -875
Out[33]= {9377.00000000000, 1.84799633554 10    , 3.94472940532 10    }

                  
Out[32]= wang-math-plot-6.png 
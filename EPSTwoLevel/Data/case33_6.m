function case33_6=case33_6()

case33_6.baseMVA=1;

case33_6.acline=[...
35	2	0.000575259	0.000297612	0	1	0	0	0	0	1	1
2	37	0.003075952	0.001566676	0	1	0	0	0	0	1	1
37	40	0.002283567	0.001162997	0	1	0	0	0	0	1	1
40	5	0.002377779	0.001211039	0	1	0	0	0	0	1	1
5	6	0.005109948	0.004411152	0	1	0	0	0	0	1	1
6	39	0.001167988	0.00386085	0	1	0	0	0	0	1	1
39	8	0.010677857	0.007706101	0	1	0	0	0	0	1	1
8	9	0.00642643	0.004617047	0	1	0	0	0	0	1	1
9	10	0.006488823	0.004617047	0	1	0	0	0	0	1	1
10	11	0.001226637	0.000405551	0	1	0	0	0	0	1	1
11	12	0.002335976	0.00077242	0	1	0	0	0	0	1	1
12	13	0.009159223	0.007206337	0	1	0	0	0	0	1	1
13	14	0.003379179	0.004447963	0	1	0	0	0	0	1	1
14	41	0.003687398	0.003281847	0	1	0	0	0	0	1	1
41	16	0.004656354	0.003400393	0	1	0	0	0	0	1	1
16	17	0.008042397	0.010737754	0	1	0	0	0	0	1	1
17	18	0.004567133	0.003581331	0	1	0	0	0	0	1	1
2	19	0.001023237	0.000976443	0	1	0	0	0	0	1	1
19	20	0.009385084	0.008456683	0	1	0	0	0	0	1	1
20	21	0.002554974	0.002984859	0	1	0	0	0	0	1	1
21	22	0.004423006	0.005848052	0	1	0	0	0	0	1	1
37	23	0.002815151	0.001923562	0	1	0	0	0	0	1	1
23	24	0.005602849	0.004424254	0	1	0	0	0	0	1	1
24	25	0.005590371	0.00437434	0	1	0	0	0	0	1	1
6	26	0.001266568	0.000645139	0	1	0	0	0	0	1	1
26	27	0.001773196	0.00090282	0	1	0	0	0	0	1	1
27	28	0.006607369	0.00582559	0	1	0	0	0	0	1	1
28	29	0.005017607	0.004371221	0	1	0	0	0	0	1	1
29	42	0.003166421	0.001612847	0	1	0	0	0	0	1	1
42	31	0.006079528	0.006008401	0	1	0	0	0	0	1	1
31	32	0.001937288	0.002257986	0	1	0	0	0	0	1	1
32	33	0.002127585	0.003308052	0	1	0	0	0	0	1	1];

%Supply node is a swing bus!
case33_6.bus=[...
24	1	0	0	0	0.42	0.2	    0	0	3	0	0	12.66	1.2	0.7
13	1	0	0	0	0.06	0.035	0	0	3	0	0	12.66	1.2	0.7
35	1	0	0	0	0	       0	0	0	3	0	0	12.66	1.2	0.7
40	1	0	0	0	0.12	0.08	0	0	3	0	0	12.66	1.2	0.7
32	1	0	0	0	0.21	0.1	    0	0	3	0	0	12.66	1.2	0.7
5	1	0	0	0	0.06	0.03	0	0	3	0	0	12.66	1.2	0.7
22	1	0	0	0	0.09	0.04	0	0	3	0	0	12.66	1.2	0.7
6	1	0	0	0	0.06	0.02	0	0	3	0	0	12.66	1.2	0.7
17	1	0	0	0	0.06	0.02	0	0	3	0	0	12.66	1.2	0.7
39	1	0	0	0	0.2	     0.1	0	0	3	0	0	12.66	1.2	0.7
19	1	0	0	0	0.09	0.04	0	0	3	0	0	12.66	1.2	0.7
8	1	0	0	0	0.2	     0.1	0	0	3	0	0	12.66	1.2	0.7
29	1	0	0	0	0.12	0.07	0	0	3	0	0	12.66	1.2	0.7
9	1	0	0	0	0.06	0.02	0	0	3	0	0	12.66	1.2	0.7
11	1	0	0	0	0.045	0.03	0	0	3	0	0	12.66	1.2	0.7
2	1	0	0	0	0.1	    0.06	0	0	3	0	0	12.66	1.2	0.7
12	1	0	0	0	0.06	0.035	0	0	3	0	0	12.66	1.2	0.7
14	1	0	0	0	0.12	0.08	0	0	3	0	0	12.66	1.2	0.7
41	1	0	0	0	0.06	0.01	0	0	3	0	0	12.66	1.2	0.7
26	1	0	0	0	0.06	0.025	0	0	3	0	0	12.66	1.2	0.7
16	1	0	0	0	0.06	0.02	0	0	3	0	0	12.66	1.2	0.7
18	1	0	0	0	0.09	0.04	0 -0.1	3	0	0	12.66	1.2	0.7
20	1	0	0	0	0.09	0.04	0	0	3	0	0	12.66	1.2	0.7
37	1	0	0	0	0.09	0.04	0	0	3	0	0	12.66	1.2	0.7
21	1	0	0	0	0.09	0.04	0	0	3	0	0	12.66	1.2	0.7
23	1	0	0	0	0.09	0.05	0	0	3	0	0	12.66	1.2	0.7
25	1	0	0	0	0.42	0.2	    0	0	3	0	0	12.66	1.2	0.7
27	1	0	0	0	0.06	0.025	0	0	3	0	0	12.66	1.2	0.7
10	1	0	0	0	0.06	0.02	0	0	3	0	0	12.66	1.2	0.7
28	1	0	0	0	0.06	0.02	0	0	3	0	0	12.66	1.2	0.7
42	1	0	0	0	0.2	     0.6	0	0	3	0	0	12.66	1.2	0.7
31	1	0	0	0	0.15	0.07	0	0	3	0	0	12.66	1.2	0.7
33	1	0	0	0	0.06	0.04	0	0	1	0	0	12.66	1.2	0.7];
				

classdef Decoder
    methods(Static)
         function tps = Exctract_TPS(frame)

            TPS_K = [34 50 209 346 413 569 595 688 790 901 1073 1219 1262 1286 1469 1594 1687];
            TPS_symb = frame(:, TPS_K+1);
            tps = -1*ones(67,1);
            for j=2:68

                if mean(abs(TPS_symb(j,:) - TPS_symb(j-1,:))) < 1e-3
                    tps(j-1) = 0;
                else 
                    tps(j-1) = 1;

                end
            end

    
         end

         function bits = QAM2Bi(ak)

             ak = round(ak*sqrt(42));
             bits = zeros(length(ak),6);

             ak_r = real(ak)+8;
             ak_i = imag(ak)+8;

             map = NaN(1, 15);
             map([1, 3, 5, 7, 9, 11, 13, 15]) = [4, 5, 7, 6, 2, 3, 1, 0];

             sym_r = map(ak_r);
             sym_i = map(ak_i);

             bits(:,[1,3,5]) = de2bi(sym_r, 'left-msb');
             bits(:,[2,4,6]) = de2bi(sym_i, 'left-msb');

         end


         function qamsymb = Serial_Data(data)
             qamsymb = data.';
             qamsymb = qamsymb(qamsymb ~= 0);

         end


         function data = Exctract_Pilots(symbols, Kmax)
             k_continual_pilot = [0 48 54 87 141 156 192 201 255 279 282 333 432 450 483 525 531 618 636 714 759 765 780 804 873 888 918 939 942 969 984 1050 1101 1107 1110 1137 1140 1146 1206 1269 1323 1377 1491 1683 1704];
             TPS_K = [34 50 209 346 413 569 595 688 790 901 1073 1219 1262 1286 1469 1594 1687];

             num_symbols = size(symbols,2);
             for i = 1:num_symbols
                 symbols(i, k_continual_pilot + 1) = 0;
                 symbols(i, TPS_K + 1) = 0;

                 k_scatter_pilot = 3*mod(i-1,4):12:Kmax;
                 symbols(i, k_scatter_pilot + 1) = 0;
             end
             data = symbols;

         end

         function idx_sync = Find_Sync(TPS)
             idx_sync = 0;
             TPS = TPS(:).';
             sync_tps = [[0,0,1,1,0,1,0,1,1,1,1,0,1,1,1,0];
                 [1,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1]];

             sync_1 = strfind(TPS, sync_tps(1,:));
             sync_2 = strfind(TPS, sync_tps(2,:));

             if isequal(size(sync_1), [0,0])
                 if isequal(size(sync_2), [0,0])
                     return
                 else
                     idx_sync = sync_2;
                 end
             else
                 idx_sync = sync_1;
             end


         end


         function deintrlv_symb = Symb_Deintrlv(intrl_data)

            N_max = 1512; % 2k mode
            deintrlv_symb = NaN;
            if mod(size(intrl_data,2), N_max) ~= 0
                return;
            end

            H = [1	1025	513	1029	33	1281	9	1153	2	1027	21	1121	769	1037	161	1282	11	1173	98	25	1217	514	1031	53	1377	777	1165	162	1284	31	1269	866	185	1474	524	1171	86	785	1101	673	1286	43	1429	106	26	1219	534	1127	821	1389	937	1422	172	1432	128	890	698	1480	576	862	946	1360	703	844	210	531	1105	577	37	1313	265	1161	130	1028	23	1141	865	173	1442	268	1183	246	795	1241	706	51	1365	617	46	1443	286	1259	918	1136	951	1402	971	200	371	749	304	1471	510	924	1244	728	851	709	295	1341	489	144	1208	376	1017	688	1460	352	1010	699	1490	588	82	529	1093	545	1285	41	1409	10	1155	22	1123	789	1133	929	1294	171	1430	108	122	538	1223	566	1383	829	938	1424	192	896	986	692	1364	607	838	433	1354	651	1170	68	113	525	1189	290	1291	157	1250	772	1051	213	775	1081	449	135	1078	355	237	272	1215	502	923	1242	708	83	613	301	1449	394	1164	152	1144	887	973	424	1344	511	1008	476	756	603	582	305	1353	649	1158	36	1303	125	782	1199	438	1388	927	1274	964	211	615	333	390	1068	407	1150	995	207	360	505	656	1204	344	1009	687	1458	332	242	539	1233	578	49	1345	521	1157	34	1283	29	1249	770	1039	181	1378	779	1177	194	19	1109	609	45	1441	266	1163	150	1124	791	1145	961	167	1334	363	238	284	1247	758	827	1497	714	52	1367	637	814	1455	446	928	1276	984	979	711	327	485	399	1214	484	251	624	346	658	1096	567	1393	841	166	1316	287	1277	994	187	1494	620	90	530	1095	565	1381	809	1421	170	1412	32	1271	886	953	1486	684	1428	96	882	697	1478	556	1427	94	786	1103	693	1382	811	1433	202	20	1111	629	813	1453	426	1420	160	1272	888	985	680	1332	351	998	443	1502	748	92	626	569	1477	554	1415	62	1507	798	1263	950	1392	959	972	212	627	589	294	1323	413	1258	900	1052	215	871	461	392	1088	503	1007	464	500	763	592	338	657	1094	547	1297	73	6	1059	277	1129	897	1038	163	1302	107	110	282	1227	662	1128	823	1401	969	168	1336	383	1006	444	1504	768	860	722	563	1361	585	38	1315	285	1257	898	1040	183	1398	875	206	276	1119	757	815	1465	458	148	1112	631	845	422	1324	415	1278	996	219	616	345	646	1064	311	1405	1001	176	1464	384	1018	700	1492	608	850	689	1350	555	1425	74	18	1091	533	1125	801	1293	169	1410	12	1175	118	793	1229	674	1288	63	874	186	1476	544	1267	854	945	1358	683	1426	76	114	537	1221	546	1287	61	1505	778	1167	182	1380	799	1273	962	179	1366	619	78	274	1099	661	1126	803	1305	201	8	1079	373	909	1198	420	1312	255	880	474	660	1108	599	837	421	1322	395	1182	228	123	622	314	1483	670	1256	824	1403	989	936	1340	479	1000	475	744	347	742	315	1501	746	60	1495	638	826	1487	702	1512	832	990	948	1372	735	840	465	647	1074	323	229	271	1213	482	155	1238	612	89	518	1063	309	1385	905	1166	164	1304	127	878	442	1484	672	1268	856	977	679	1330	331	230	283	1245	738	59	1493	618	58	1475	542	1255	822	1391	957	940	1436	224	884	729	552	1331	349	902	1072	439	1406	1003	208	372	761	560	1459	350	914	1104	695	1394	843	198	275	1117	737	47	1461	362	154	1220	536	1139	853	933	1326	427	1438	236	124	634	570	1479	574	1511	830	958	960	992	980	723	583	325	389	1066	387	1054	227	111	366	410	1228	664	1140	855	965	423	1342	491	240	380	762	572	1491	606	818	1359	701	1510	812	1435	222	788	1115	725	807	1337	457	136	1080	375	1005	432	1472	512	1020	732	596	593	549	1317	297	1417	138	1156	24	1143	885	941	1454	428	1440	256	892	730	564	1363	605	806	1327	445	908	1180	216	883	717	296	1343	509	912	1212	472	1011	719	328	497	655	1202	324	241	527	1201	322	145	1090	515	1041	65	5	1057	257	1033	129	1026	3	1045	97	13	1185	258	1035	149	1122	771	1049	193	7	1077	353	141	1186	260	1055	245	783	1209	450	147	1110	611	77	262	1067	405	1130	899	1050	195	103	365	398	1196	408	1152	1015	975	456	499	751	336	498	667	1234	580	81	517	1061	289	1289	137	1154	4	1047	117	781	1197	418	1292	159	1270	868	217	520	1075	341	901	1070	419	1310	235	112	378	666	1224	568	1395	861	934	1328	447	1004	220	628	601	550	1319	317	906	1168	184	1400	895	974	436	1376	767	848	466	659	1106	579	69	261	1065	385	1034	131	1046	99	109	270	1195	406	1132	919	1146	963	199	359	493	400	1216	504	1019	720	340	753	559	1457	330	146	1092	535	1137	833	165	1314	267	1181	226	27	1237	610	57	1473	522	1159	54	1379	797	1261	930	1296	191	876	218	532	1107	597	805	1325	425	1418	140	1176	120	889	686	1448	320	1022	956	1500	736	852	721	551	1329	329	134	1060	279	1149	993	175	1462	364	250	540	1235	598	817	1357	681	1414	44	1431	126	794	1231	694	1384	831	970	180	1368	639	846	434	1356	671	1266	836	209	519	1073	321	133	1058	259	1053	225	15	1205	354	153	1218	516	1043	85	773	1069	417	1290	139	1174	100	121	526	1191	310	1387	925	1262	932	1308	223	872	473	648	1076	343	997	431	1470	492	252	636	602	562	1351	573	1509	810	1423	190	1508	800	1275	982	947	1370	715	72	369	653	1190	292	1311	253	784	1211	470	915	1114	707	71	357	397	1194	388	1056	247	879	462	404	1120	759	847	454	403	1118	739	79	358	409	1226	644	1044	87	869	429	1450	396	1184	248	891	718	308	1375	765	816	1467	478	916	1116	727	839	453	391	1086	483	239	368	506	668	1236	600	849	677	1318	299	1437	234	28	1239	630	825	1485	682	1416	64	894	954	1488	704	864	978	691	1362	587	70	273	1097	641	1030	35	1301	105	14	1187	278	1131	917	1134	931	1306	203	104	377	654	1192	312	1407	1021	944	1468	480	1012	731	584	337	645	1062	291	1309	233	16	1207	374	921	1230	676	1300	95	870	441	1482	652	1172	88	881	685	1446	300	1439	254	796	1243	726	819	1369	713	40	1335	381	910	1200	440	1408	1023	976	468	755	591	326	401	1098	643	1042	67	101	269	1193	386	1036	151	1142	867	205	264	1087	501	911	1210	452	243	623	334	402	1100	663	1138	835	197	263	1085	481	143	1206	356	249	528	1203	342	913	1102	675	1298	75	102	281	1225	642	1032	55	1397	873	174	1444	288	1279	1014	955	1498	716	84	625	557	1445	298	1419	158	1252	792	1147	981	935	1338	459	232	379	750	316	1503	766	828	1499	734	820	1371	733	808	1339	477	904	1084	471	999	463	488	507	752	348	754	571	1489	586	50	1347	541	1253	802	1295	189	1506	780	1179	214	787	1113	705	39	1333	361	142	1188	280	1151	1013	943	1466	460	244	635	590	306	1355	669	1254	804	1307	221	776	1083	469	903	1082	451	231	367	494	412	1248	760	859	710	307	1373	745	48	1463	382	922	1232	696	1396	863	966	435	1374	747	80	370	665	1222	548	1299	93	774	1071	437	1386	907	1178	196	115	621	302	1451	414	1260	920	1148	983	967	455	487	495	496	508	764	604	594	561	1349	553	1413	42	1411	30	1251	790	1135	949	1390	939	1434	204	116	633	558	1447	318	926	1264	952	1404	991	968	467	743	335	486	411	1246	740	91	614	313	1481	650	1160	56	1399	893	942	1456	448	1024	988	724	595	581	293	1321	393	1162	132	1048	119	877	430	1452	416	1280	1016	987	712	339	741	303	1469	490	156	1240	632	857	678	1320	319	1002	188	1496	640	858	690	1352	575	842	178	1348	543	1265	834	177	1346	523	1169	66	17	1089];
            m = size(intrl_data,2)/N_max;
            deintrlv_symb = 3*zeros(size(intrl_data));
            c = 0;
            for k = 0:N_max:m*N_max-1
                idx = (k+1):(k+N_max);
                if mod(c,2) == 0
                    deintrlv_symb(:,idx) = intrl_data(:, H + c*N_max);
                else
                     deintrlv_symb(:,H + c*N_max) = intrl_data(:,idx);
                end
                c= c+1;   
            end

         end
         

         function deintrlv_bits = Bits_Deintrlv(intrl_bits)
             deintrlv_bits = NaN;
             if mod(size(intrl_bits,2), 126) ~= 0
                 return;
             end


             H_bit = zeros(6, 126);
             % Funciones de permutacion
             H_bit(1, :) = (0:125) + 1;
             H_bit(2, :) = mod((0:125) + 63, 126) + 1;
             H_bit(3, :) = mod((0:125) + 105, 126)+ 1;
             H_bit(4, :) = mod((0:125) + 42, 126)+ 1;
             H_bit(5, :) = mod((0:125) + 21, 126)+ 1;
             H_bit(6, :) = mod((0:125) + 84, 126)+ 1;

             m = size(intrl_bits,2)/126;

             idx = zeros(size(H_bit,1), m*size(H_bit,2));

             for i = 0:m-1
                 idx(:, 126*i+1: 126*i+126) =  H_bit + 126*i;
             end
             deintrlv_bits = zeros(size(idx));
             deintrlv_bits(1,idx(1,:)) = intrl_bits(1,:);
             deintrlv_bits(2,idx(2,:)) = intrl_bits(2,:);
             deintrlv_bits(3,idx(3,:)) = intrl_bits(3,:);
             deintrlv_bits(4,idx(4,:)) = intrl_bits(4,:);
             deintrlv_bits(5,idx(5,:)) = intrl_bits(5,:);
             deintrlv_bits(6,idx(6,:)) = intrl_bits(6,:);

             deintrlv_bits = deintrlv_bits(:);
         end

            

    end
end
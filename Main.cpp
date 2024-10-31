#include "Param.h"

#include <fstream>
#include <iostream>
#include <stdio.h>
#include "Simulation.h"
#include <time.h>
#include <string>

using namespace std;


int main()
{
	int lineRateType;
	ofstream outfile;
	Simulation simulation;
	Database dataBase;

	outfile.open("NID_B_ILP.xls", ios::trunc);//ios::out
	if (!outfile)
	{
		cerr << "存取数据的文件打不开" << endl;
	}
	else
	{
		outfile << "SFCNum \t" << "TotalRev\t" << "AveTotalCost \t" << "Aveshouyi\t" << "AveTotalSlots \t" << "AveTotalTrans \t" << "AveTotalRegen \t\t"
			<< "totalS\t" << "BP\t" << "TIME\t" << endl;
		outfile.close();
	}

	outfile.open("BP-N-kun-a=0.8-ILP.xls", ios::trunc);
	if (!outfile)
	{
		cerr << "存取数据的文件打不开" << endl;
	}
	else
	{
		outfile.close();
	}

	outfile.open("SP-NID-kun.xls", ios::trunc);
	if (!outfile)
	{
		cerr << "存取数据的文件打不开" << endl;
	}
	else
	{
		outfile.close();
	}

	//估计是设置功能链的长度
	for (virMaxNode = 3; virMaxNode < 4; virMaxNode++)
	{
		for (int ithTimeSFC = 0; ithTimeSFC < runTime; ithTimeSFC++)
		{
			simulation.minMatrix();
			SFCnum = totalSFC;

			minNodesWeight = INFValue;
			conNumBreakEvent = 0;
			rejectConnectionRequestNum = 0;
			acceptNum = 0;
			totalSpectrumResource = 0;

			totalGBbs = 0;
			totalEnergy = 0;
			totalSpectrum = 0;
			TRANum = 0;
			REGNum = 0;
			WASTEBAND = 0;
			printer = 1;
			rejectConnectionRequestNum = 0;
			acceptNum = 0;
			totalSpectrumResource = 0;
			totalSpectrum = 0;
			totalGBbs = 0;
			totalSpectrum = 0;
			/*totalSpectrum1 = 0;*/

			TRANum = 0;
			REGNum = 0;

			TRANum400_1500 = 0;
			TRANum100_1700 = 0;
			TRANum40_1800 = 0;

			REGNum400_1500 = 0;
			REGNum100_1700 = 0;
			REGNum40_1800 = 0;
			//WASTEBAND = 0;

			minEnergy = INFValue;
			int ithSFC = 0;
			srand(ithTimeSFC);
			//srand(1);
			while (ithSFC < SFCnum)//SFCnum=totalSFC，业务数
			{
				simulation.m_generateEvent. y(ithSFC);
				ithSFC++;
			}
			cout << ithTimeSFC << "kun" << endl;
			simulation.mapping(ithTimeSFC, SFCnum);	//估计是映射函数
			simulation.m_database.resourceUtilizationRatio();
		}
		simulation.m_database.statisticResult();

		if ((int)simulation.m_database.allResults.size() != 0)
		{
			simulation.m_database.allResults.clear();
		}
	}

	system("pause");
	return 0;
}

/*
int main()
{
	//int lineRateType;
	Database database;
	Simulation simulation;
	//Suurballe suurballe;

	ofstream outfile;
	//ofstream outfile;
	outfile.open("Result_MinEnergy_Suurballe_NSFNET_EnergyMatrix.xls", ios::trunc);
	if(!outfile)
	{
		cerr << "存取数据的文件打不开" << endl;
	}
	else
	{
		outfile<<"K_SP\t"<<"con_Res_Num\t"<<"BP\t"<<"slotNum\t"<<"aveEnergyGBps\t"
			<<"aveTran\t"<<"aveRegen\t"<<"aveWasteBand\t"<<endl;
		outfile.close();
	}

	//database.graph.readTopology(MaxNodeNum,"NSF_KM.txt");
	database.graph.readTopology(MaxNodeNum,"USNET_KM.txt");
	//database.graph.readTopology(MaxNodeNum,"75NodeWeightedTopology.txt");
	//srand(2);
	//database.graph.readTopology(MaxNodeNum,"6nodeAllConnection.txt");
	//database.graph.readTopology(MaxNodeNum,"NSFtopology.txt");
	//m_database.graph.readTopology(MaxNodeNum,"6node-topology.txt");
	//database.graph.failProbability();
	int allStatus = 0;
	int guardBand = 0;
	while(guardBand < GBNum)
	{
		GB = guardBand;
		outfile.open("Result_MinEnergy_Suurballe_NSFNET_EnergyMatrix.xls", ios::app);
		if(!outfile)
		{
			cerr << "存取数据的文件打不开" << endl;
		}
		else
		{
			outfile<<"GB = "<<guardBand<<endl;
			outfile.close();
		}


		eachCRNum = EACHCRNum;
		acceptNumID = 0;
		int ith = 0;
		for (int k = 0; k < 10; k++)
		{
			conNumBreakEvent = 0;
			rejectConnectionRequestNum = 0;
			acceptNum = 0;

			totalSpectrumResource = 0;
			totalGBbs = 0;
			totalEnergy = 0;
			totalSpectrum = 0;
			TRANum = 0;
			REGNum = 0;

			srand(k);
			bool val = simulation.getTraffic();
			if (val==false)
			{
				continue;
			}
			else if(val == true)
			{
				cout<<"k:  "<<k<<endl;
				ith++;
				//database.resourceUtilizationRatio();
				//database.statisticResult();
			}
			if (ith == runTime)
			{
				break;
			}
			else if (k==100)
			{
				allStatus = 1;
				//cout<<"There are no enough spectrum for all CRs"
			}
		}

		if (allStatus == 1)
		{
			cout<<"There are no enough spectrum for all CRs"<<endl;
			break;
		}
		int j = 0;
		K_SP = 3;
		ithTime = 0;
		while (j < K_SP_Time)
		{
			start_ksp = j;
			//cout<<"j  "<<j<<endl;

			ConnectionRequestNum = CRNUM;

			//conNumBreakEvent = 0;
			for(int i = 0; i< figNodeNum; i++)
			{
				start_conNum = i;

				outfile.open("Result_MinEnergy_Suurballe_NSFNET_EnergyMatrix.xls", ios::app);
				if(!outfile)
				{
					cerr << "存取数据的文件打不开" << endl;
				}
				else
				{
					outfile<<K_SP<<"\t"<<ConnectionRequestNum<<"\t";
					outfile.close();
				}

				ithTime = i;
				conNumBreakEvent = 0;
				rejectConnectionRequestNum = 0;
				acceptNum = 0;
				totalSpectrumResource = 0;

				totalGBbs = 0;
				totalEnergy = 0;
				totalSpectrum = 0;
				TRANum = 0;
				REGNum = 0;
				WASTEBAND = 0;
				printer = 1;

				for (virMaxNode = 3; virMaxNode< 7; virMaxNode++)
				{
					minNodesWeight = INFValue;

					totalGBbs = 0;
					totalEnergy = 0;
					PRIEnergy = 0;
					PROEnergy = 0;
					TRANEnergy = 0;
					RENEnergy = 0;

					totalSpectrum = 0;
					PRISpectrum = 0;
					PROSpectrum = 0;

					TRANum = 0;
					REGNum = 0;
					PRIRegenNum = 0;
					PRORegenNum = 0;

					totalHopNum = 0;
					PRIHopNum = 0;
					PROHopNum = 0;
					simulation.minMatrix();
					for (int k=0; k<virTopoTime;k++)
					{
						srand(k);
						minEnergy = INFValue;
						simulation.m_graph.ERtopology();
						simulation.mapping();
						//simulation.run();
					}
					simulation.statisticResult();
				}

				for (int k = 0; k< runTime; k++)
				{
					simulation.run(k);
				}
				database.statisticResult();

				ConnectionRequestNum = ConnectionRequestNum + GAP;
			}
			K_SP++;
			j++;
		}
		while(!Event_pq_tem.empty())
		{
			Event_pq_tem.top();
			Event_pq_tem.pop();
		}
		guardBand++;
	}
	return 0;
}*/
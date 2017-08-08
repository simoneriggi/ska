// **************************************************************************
// * License and Disclaimer                                                 *
// *                                                                        *
// * Copyright 2016 Simone Riggi																			      *
// *																																	      *
// * This file is part of SKA DSH.LMC 																		  *
// * SKA DSH.LMC is free software: you can redistribute it and/or modify it *
// * under the terms of the GNU General Public License as published by      *
// * the Free Software Foundation, either * version 3 of the License,       * 
// * or (at your option) any later version.                                 * 
// * Caesar is distributed in the hope that it will be useful, but 			    *
// * WITHOUT ANY WARRANTY; without even the implied warranty of             * 
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
// * See the GNU General Public License for more details. You should        * 
// * have received a copy of the GNU General Public License along with      * 
// * Caesar. If not, see http://www.gnu.org/licenses/.                      *
// **************************************************************************

/**
* @file Test_StatsUtils.cc
* @class Test_StatsUtils
* @brief Class for testing StatsUtils class
*
* Test class for StatsUtils class
* @author S. Riggi
* @date 08/08/2017
*/


#include <Test_StatsUtils.h>
#include <StatsUtils.h>

#include "gtest/gtest.h"
#include "gmock/gmock-matchers.h"

#include <iostream>


namespace Caesar {

//Set static vars
std::vector<double> StatsUtilsTest::m_data;
double StatsUtilsTest::m_numeric_tol= 1.e-6;

//===========================================
//==           DEFINE TESTS              ====
//===========================================

TEST_F(StatsUtilsTest, TestMedianFast)
{
	
	//Expected value
	//NB: Using data generated in astropy example (ref. http://docs.astropy.org/en/stable/stats/robust.html#sigma-clipping)
	double median_exp= 0.032658645;//computed using python numpy
	//double median_exp= 0.03265865;//computed using R

	//Compute median
	bool useParallelVersion= false;
	double median= Caesar::StatsUtils::GetMedianFast(m_data,useParallelVersion);
	double diff= fabs(median-median_exp);

	//Compute median with parallel algorithms
	useParallelVersion= true;
	double median_parallel= Caesar::StatsUtils::GetMedianFast(m_data,useParallelVersion);
	double diff_parallel= fabs(median_parallel-median_exp);

	//Print text
	std::stringstream ss;
	ss<<"StatsUtilsTest::TestMedianFast(): INFO: median_exp="<<median_exp<<", median="<<median<<" (diff="<<diff<<"), median_parallel="<<median_parallel<<" (diff="<<diff_parallel<<")";
	std::cerr<<ss.str()<<std::endl;

	//Check if median computation is correct within numerical tolerance
	ASSERT_LT(diff,m_numeric_tol);
	ASSERT_LT(diff_parallel,m_numeric_tol);
	

}//close TestMedianFast()



TEST_F(StatsUtilsTest, TestMADFast)
{
	double mad_exp= 0.30728245008874866;//computed using python mad_std in astropy
	//double mad_exp= 0.307282;//computed using R
	
	//Compute mad
	bool useParallelVersion= false;
	double median= Caesar::StatsUtils::GetMedianFast(m_data,useParallelVersion);
	double mad= Caesar::StatsUtils::GetMADFast(m_data,median,useParallelVersion);
	mad*= 1.4826;
	double diff= fabs(mad-mad_exp);

	//Compute mad with parallel algorithms
	useParallelVersion= true;
	double median_parallel= Caesar::StatsUtils::GetMedianFast(m_data,useParallelVersion);
	double mad_parallel= Caesar::StatsUtils::GetMADFast(m_data,median_parallel,useParallelVersion);
	mad_parallel*= 1.4826;
	double diff_parallel= fabs(mad_parallel-mad_exp);

	//Print text
	std::stringstream ss;
	ss<<"StatsUtilsTest::TestMADFast(): INFO: mad_exp="<<mad_exp<<", mad="<<mad<<" (diff="<<diff<<"), mad_parallel="<<mad_parallel<<" (diff="<<diff_parallel<<")";
	std::cerr<<ss.str()<<std::endl;

	//Check if mad computation is correct within numerical tolerance
	ASSERT_LT(diff,m_numeric_tol);
	ASSERT_LT(diff_parallel,m_numeric_tol);
	

}//close TestMADFast()

TEST_F(StatsUtilsTest, TestMeanRMS)
{
	double mean_exp= 0.86586417685994999;//computed using numpy mean
	double rms_exp= 3.2996406358899497;//computed using numpy std(data,ddof=2)

	//Compute mean & rms
	double mean= 0;
	double rms= 0;
	Caesar::StatsUtils::ComputeMeanAndRMS(mean,rms,m_data);
	double diff_mean= fabs(mean-mean_exp);
	double diff_rms= fabs(rms-rms_exp);
	
	//Compute mean & rms with parallel versions
	double mean_parallel= 0;
	double rms_parallel= 0;
	Caesar::StatsUtils::ComputeMeanAndRMS(mean_parallel,rms_parallel,m_data);
	double diff_mean_parallel= fabs(mean_parallel-mean_exp);
	double diff_rms_parallel= fabs(rms_parallel-rms_exp);

	//Print text
	std::stringstream ss;
	ss<<"StatsUtilsTest::TestMeanRMS(): INFO: mean_exp="<<mean_exp<<", mean="<<mean<<" (diff="<<diff_mean<<"), mean_parallel="<<mean_parallel<<" (diff="<<diff_mean_parallel<<"), rms_exp="<<rms_exp<<", rms="<<rms<<" (diff="<<diff_rms<<"), rms_parallel="<<rms_parallel<<" (diff="<<diff_rms_parallel<<")";
	std::cerr<<ss.str()<<std::endl;

	//Check if mean/rms computation is correct within numerical tolerance
	ASSERT_LT(diff_mean,m_numeric_tol);
	ASSERT_LT(diff_mean_parallel,m_numeric_tol);
	ASSERT_LT(diff_rms,m_numeric_tol);
	ASSERT_LT(diff_rms_parallel,m_numeric_tol);

}//close TestMeanRMS()

TEST_F(StatsUtilsTest, TestClippedEstimators)
{

	//Expected values @ sigmaclip=3, niter= 10 (computed using python astropy.stats.sigma_clipped_stats()
	double mean_clipped_exp= -0.0020337793885611471;
	double median_clipped_exp= -0.023632809000000001;
	double rms_clipped_exp= 0.19585230166658546;

	//Expected value @ sigmaclip=3, niter= 1
	double mean_clipped_exp_0= 0.44795528741963725;
	double median_clipped_exp_0= 0.018870317899999999;
	double rms_clipped_exp_0= 2.467320281271788;

	//Compute median/mean/rms
	bool useParallelVersion= false;
	double mean= 0;
	double rms= 0;
	Caesar::StatsUtils::ComputeMeanAndRMS(mean,rms,m_data);
	double median= Caesar::StatsUtils::GetMedianFast(m_data,useParallelVersion);

	//Compute clipped estimators
	ClippedStats<double> stats;
	double sigmaclip= 3;
	int niter= 10;
	double tol= 1.e-6;
	Caesar::StatsUtils::GetClippedEstimators(stats,m_data,median,mean,rms,sigmaclip,niter,tol,useParallelVersion);
	
	double mean_clipped= stats.mean;
	double median_clipped= stats.median;
	double rms_clipped= stats.stddev;
	
	double diff_mean= fabs(mean_clipped-mean_clipped_exp);
	double diff_median= fabs(median_clipped-median_clipped_exp);
	double diff_rms= fabs(rms_clipped-rms_clipped_exp);

	//Print text
	std::stringstream ss;
	ss<<"StatsUtilsTest::TestClippedEstimators(): INFO: mean_exp="<<mean_clipped_exp<<", mean="<<mean_clipped<<" (diff="<<diff_mean<<"), rms_clipped_exp="<<rms_clipped_exp<<", rms="<<rms_clipped<<" (diff="<<diff_rms<<"), median_exp="<<median_clipped_exp<<", median="<<median_clipped<<" (diff="<<diff_median<<")";
	std::cerr<<ss.str()<<std::endl;

	//Check if computation is correct within numerical limits
	ASSERT_LT(diff_mean,m_numeric_tol);
	ASSERT_LT(diff_median,m_numeric_tol);
	ASSERT_LT(diff_rms,m_numeric_tol);

}//close TestClippedEstimators()


}//close namespace



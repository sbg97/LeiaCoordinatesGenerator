#define PI 3.141592654
#define SECONDS_PER_POINT 0.001
#define FPS 15
#include <algorithm>
#include <iostream>
#include <math.h>
#include <vector>
#include "qcustomplot/qcustomplot.h"

using namespace std;

// takes the times and points from input_times and input_coord
// interpolates them using Hermite
// and puts them in output_coord
void perform_interpolation(const vector<double> *input_times, const vector<double> *input_coord, vector<double> *output_coord)
{

	output_coord->clear();

	double current_time = input_times->at(0);
	while (current_time <= input_times->back())
	{
		// find p0, p1, m0, and m1

		// find the index of the previous point
		size_t k_0 = 0;
		bool equal_done = false;

		double interpolated_point;

		for (size_t i = 0; i < input_times->size(); i++)
		{
			if (input_times->at(i) == current_time)
			{
				// just put the point in and be done with it
				interpolated_point = input_coord->at(i);
				equal_done = true;
				break;
			}
			if (input_times->at(i) > current_time)
			{
				k_0 = i - 1;
				break;
			}
		}
		if (!equal_done)
		{
			// find the index of the next point
			size_t k_1 = k_0 + 1;

			// find the slope of the previous point
			double m_0;
			// edge case when k_0 == 0
			if (k_0 == 0)
			{
				m_0 = (input_coord->at(k_1) - input_coord->at(k_0)) / (input_times->at(k_1) - input_times->at(k_0));
			}
			// k_0>0, so k_0-1 exists
			else
			{
				m_0 = 0.5 * ((input_coord->at(k_1) - input_coord->at(k_0)) / (input_times->at(k_1) - input_times->at(k_0)) + (input_coord->at(k_0) - input_coord->at(k_0 - 1)) / (input_times->at(k_0) - input_times->at(k_0 - 1)));
			}

			// find the slope of the next point
			double m_1;
			// edge case when k_1 == input_times->size() - 1
			if (k_1 == input_times->size() - 1)
			{
				m_1 = (input_coord->at(k_1) - input_coord->at(k_0)) / (input_times->at(k_1) - input_times->at(k_0));
			}
			// k_1<input_times->size() - 1, so k_1+1 exists
			else
			{
				m_1 = 0.5 * ((input_coord->at(k_1 + 1) - input_coord->at(k_1)) / (input_times->at(k_1 + 1) - input_times->at(k_1)) + (input_coord->at(k_1) - input_coord->at(k_0)) / (input_times->at(k_1) - input_times->at(k_0)));
			}

			// find the distance we need to interpolate as a double between 0 and 1, 0 being point 0, and 1 being point 1
			double t = (current_time - input_times->at(k_0)) / (input_times->at(k_1) - input_times->at(k_0));

			// use t to get h00, h10, h01, and h11
			double h00 = 2 * t * t * t - 3 * t * t + 1;
			double h10 = t * t * t - 2 * t * t + t;
			double h01 = -2 * t * t * t + 3 * t * t;
			double h11 = t * t * t - t * t;

			// try renormalizing m_0 and m_1
			m_0 *= (input_times->at(k_1) - input_times->at(k_0));
			m_1 *= (input_times->at(k_1) - input_times->at(k_0));

			interpolated_point = h00 * input_coord->at(k_0) + h10 * m_0 + h01 * input_coord->at(k_1) + h11 * m_1;
		}

		// put in the interpolated coordinate
		output_coord->push_back(interpolated_point);

		// continue with the next time step
		current_time += SECONDS_PER_POINT;
	}
}

void clamp0to255(vector<double> *output_coord)
{
	for (size_t i = 0; i < output_coord->size(); i++)
	{
		if (output_coord->at(i) > 255)
		{
			output_coord->at(i) = 255;
		}
		if (output_coord->at(i) < 0)
		{
			output_coord->at(i) = 0;
		}
	}
}

int main()
{
	vector<double> xcoords;
	vector<double> ycoords;
	vector<double> zcoords;
	vector<double> camcoords;
	vector<double> rvals;
	vector<double> gvals;
	vector<double> bvals;
	vector<double> time_coords;

	double explosionxcoords[] = {0, 1, -5, -2, -8, -3, -7, -1};
	double explosionycoords[] = {0, 4, 8, 1, 0, -2, -7, -3};
	double explosion2xcoords[] = {0, 1, -5, -2, -8, -3, -7, -1, -1, 2};
	double explosion2ycoords[] = {0, 4, 8, 1, 0, -2, -7, -3, -8, -4};

	//Shoot 2 green disruptors from (206, 85, 70) to (64, 47, 70) //z used to be 55
	//Should take 1.8 seconds @ 3.5ms/point = 514 points
	{
		double start_x = 250;
		double start_y = 15;
		double start_z = 70;
		double end_x = 125;
		double end_y = 109;
		double end_z = 70;
		double disruptor_percentage = 0.40;
		double total_duration = 0.249;
		double start_time = 0.0;
		//first shot
		for (double i = start_time; i < total_duration; i += 1.0 / FPS)
		{

			double interp;

			time_coords.push_back(i);

			// insert close point
			interp = max((1.0 + disruptor_percentage) / total_duration * (time_coords.back() - start_time) - disruptor_percentage, 0.0);
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(255);
			bvals.push_back(0);
			camcoords.push_back(135);

			// a little time has passed between points
			time_coords.push_back(time_coords.back() + ((1.0 / FPS) / 2.0));

			// insert far point
			interp = min((1.0 + disruptor_percentage) / total_duration * (time_coords.back() - start_time), 1.0);
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(255);
			bvals.push_back(0);
			camcoords.push_back(135);
		}
		//little explosion at (64, 47, 70) with -100 deg rotation
		start_x = 125;
		start_y = 109;
		start_z = 70;
		total_duration = 0.1295;
		start_time = time_coords.back();
		double rotation = 100.0 / 180.0 * PI;
		double scale = 4.0;
		for (double i = start_time; i < start_time + total_duration; i += 1.0 / FPS)
		{
			for (int j = 0; j < 8; j++)
			{
				time_coords.push_back(time_coords.back() + 1.0 / FPS / 8.0);

				xcoords.push_back(start_x + scale * (explosionxcoords[j % 8] * cos(rotation) + explosionycoords[j % 8] * sin(rotation)));
				ycoords.push_back(start_y + scale * (explosionycoords[j % 8] * cos(rotation) - explosionxcoords[j % 8] * sin(rotation)));
				zcoords.push_back(start_z);
				rvals.push_back(255);
				gvals.push_back(0);
				bvals.push_back(0);
				camcoords.push_back(135);
			}
		}
		//travel back  (64, 47, 70) to (206, 85, 70)
		start_x = 125;
		start_y = 109;
		start_z = 70;
		end_x = 250;
		end_y = 15;
		end_z = 70;
		total_duration = 0.035;
		start_time = time_coords.back();
		for (double i = 0; i < 4; i++)
		{

			time_coords.push_back(time_coords.back() + total_duration / 4);

			double interp = (double)(time_coords.back() - start_time) / total_duration;
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(0);
			bvals.push_back(0);
			camcoords.push_back(135);
		}
		{
			//second shot
			double start_x = 250;
			double start_y = 15;
			double start_z = 70;
			double end_x = 125;
			double end_y = 109;
			double end_z = 70;
			double disruptor_percentage = 0.40;
			double total_duration = 0.249;
			double start_time = time_coords.back();
			//first shot
			for (double i = 0; i < total_duration; i += 1.0 / FPS)
			{

				double interp;

				time_coords.push_back(time_coords.back() + ((1.0 / FPS) / 2.0));

				// insert close point
				interp = max((1.0 + disruptor_percentage) / total_duration * (time_coords.back() - start_time) - disruptor_percentage, 0.0);
				xcoords.push_back(start_x * (1 - interp) + end_x * interp);
				ycoords.push_back(start_y * (1 - interp) + end_y * interp);
				zcoords.push_back(start_z * (1 - interp) + end_z * interp);
				rvals.push_back(0);
				gvals.push_back(255);
				bvals.push_back(0);
				camcoords.push_back(135);

				// a little time has passed between points
				time_coords.push_back(time_coords.back() + ((1.0 / FPS) / 2.0));

				// insert far point
				interp = min((1.0 + disruptor_percentage) / total_duration * (time_coords.back() - start_time), 1.0);
				xcoords.push_back(start_x * (1 - interp) + end_x * interp);
				ycoords.push_back(start_y * (1 - interp) + end_y * interp);
				zcoords.push_back(start_z * (1 - interp) + end_z * interp);
				rvals.push_back(0);
				gvals.push_back(255);
				bvals.push_back(0);
				camcoords.push_back(135);
			}
			//little explosion at (64, 47, 70) with -100 deg rotation
			start_x = 125;
			start_y = 109;
			start_z = 70;
			total_duration = 0.1295;
			start_time = time_coords.back();
			double rotation = 100.0 / 180.0 * PI;
			double scale = 4.0;
			for (double i = start_time; i < start_time + total_duration; i += 1.0 / FPS)
			{
				for (int j = 0; j < 8; j++)
				{
					time_coords.push_back(time_coords.back() + 1.0 / FPS / 8.0);

					xcoords.push_back(start_x + scale * (explosionxcoords[j % 8] * cos(rotation) + explosionycoords[j % 8] * sin(rotation)));
					ycoords.push_back(start_y + scale * (explosionycoords[j % 8] * cos(rotation) - explosionxcoords[j % 8] * sin(rotation)));
					zcoords.push_back(start_z);
					rvals.push_back(255);
					gvals.push_back(0);
					bvals.push_back(0);
					camcoords.push_back(135);
				}
			}
			//travel to next gun  (64, 47, 70) to (33, 62, 77)
			start_x = 125;
			start_y = 109;
			start_z = 70;
			end_x = 78;
			end_y = 171;
			end_z = 70;
			total_duration = 0.035;
			start_time = time_coords.back();
			for (double i = 0; i < 4; i++)
			{

				time_coords.push_back(time_coords.back() + total_duration / 4);

				double interp = (double)(time_coords.back() - start_time) / total_duration;
				xcoords.push_back(start_x * (1 - interp) + end_x * interp);
				ycoords.push_back(start_y * (1 - interp) + end_y * interp);
				zcoords.push_back(start_z * (1 - interp) + end_z * interp);
				rvals.push_back(0);
				gvals.push_back(0);
				bvals.push_back(0);
				camcoords.push_back(135);
			}
		}
	}

	//phaser from (33, 62, 77) to (246, 99, 77), scanning towards (255, 98, 77) //z used to be 77, 101, 134
	//Should take 4 seconds @ 3.5ms/point = 1142 points
	{
		double start_x = 78;
		double start_y = 171;
		double start_z = 77;
		double end_x = 90;
		double end_y = 230;
		double end_z = 77;
		double end_x2 = 140;
		double end_y2 = 200;
		double end_z2 = 77;
		double total_duration = 1.5;
		double start_time = time_coords.back();

		for (double i = 0; i < total_duration; i += 1.0 / FPS)
		{

			double interp_between_targets = (double)i / total_duration;

			time_coords.push_back(time_coords.back() + ((1.0 / FPS) / 2.0));

			// insert close point
			xcoords.push_back(start_x);
			ycoords.push_back(start_y);
			zcoords.push_back(start_z);
			rvals.push_back(255);
			gvals.push_back(0);
			bvals.push_back(0);
			camcoords.push_back(135);

			// a little time has passed between points
			time_coords.push_back(time_coords.back() + ((1.0 / FPS) / 2.0));

			// insert far point
			xcoords.push_back(end_x * (1 - interp_between_targets) + end_x2 * interp_between_targets);
			ycoords.push_back(end_y * (1 - interp_between_targets) + end_y2 * interp_between_targets);
			zcoords.push_back(end_z * (1 - interp_between_targets) + end_z2 * interp_between_targets);
			rvals.push_back(255);
			gvals.push_back(0);
			bvals.push_back(0);
			camcoords.push_back(135);
		}
		//travel to next gun  (77, 69, 77) to (206, 85, 70)
		start_x = 140;
		start_y = 200;
		start_z = 77;
		end_x = 150;
		end_y = 100;
		end_z = 70;
		total_duration = 0.035/2.0;
		start_time = time_coords.back();
		for (double i = 0; i < 4; i++)
		{

			time_coords.push_back(time_coords.back() + total_duration / 4);

			double interp = (double)(time_coords.back() - start_time) / total_duration;
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(0);
			bvals.push_back(0);
			camcoords.push_back(135);
		}
	//travel to next gun  (77, 69, 77) to (206, 85, 70)
		start_x = 150;
		start_y = 100;
		start_z = 77;
		end_x = 250;
		end_y = 15;
		end_z = 70;
		total_duration = 0.035/2.0;
		start_time = time_coords.back();
		for (double i = 0; i < 4; i++)
		{

			time_coords.push_back(time_coords.back() + total_duration / 4);

			double interp = (double)(time_coords.back() - start_time) / total_duration;
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(0);
			bvals.push_back(0);
			camcoords.push_back(135);
		}
	}

	////Shoot 3 green disruptors from (206, 85, 70) to (64, 47, 70) // z used to be 77
	{
		//second shot
		double start_x = 250;
		double start_y = 15;
		double start_z = 70;
		double end_x = 125;
		double end_y = 109;
		double end_z = 70;
		double disruptor_percentage = 0.40;
		double total_duration = 0.249;
		double start_time = time_coords.back();
		//first shot
		for (double i = 0; i < total_duration; i += 1.0 / FPS)
		{

			double interp;

			time_coords.push_back(time_coords.back() + ((1.0 / FPS) / 2.0));

			// insert close point
			interp = max((1.0 + disruptor_percentage) / total_duration * (time_coords.back() - start_time) - disruptor_percentage, 0.0);
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(255);
			bvals.push_back(0);
			camcoords.push_back(135);

			// a little time has passed between points
			time_coords.push_back(time_coords.back() + ((1.0 / FPS) / 2.0));

			// insert far point
			interp = min((1.0 + disruptor_percentage) / total_duration * (time_coords.back() - start_time), 1.0);
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(255);
			bvals.push_back(0);
			camcoords.push_back(135);
		}
		//little explosion at (64, 47, 70) with -100 deg rotation
		start_x = 125;
		start_y = 109;
		start_z = 70;
		total_duration = 0.1295;
		start_time = time_coords.back();
		double rotation = 100.0 / 180.0 * PI;
		double scale = 4.0;
		for (double i = start_time; i < start_time + total_duration; i += 1.0 / FPS)
		{
			for (int j = 0; j < 8; j++)
			{
				time_coords.push_back(time_coords.back() + 1.0 / FPS / 8.0);

				xcoords.push_back(start_x + scale * (explosionxcoords[j % 8] * cos(rotation) + explosionycoords[j % 8] * sin(rotation)));
				ycoords.push_back(start_y + scale * (explosionycoords[j % 8] * cos(rotation) - explosionxcoords[j % 8] * sin(rotation)));
				zcoords.push_back(start_z);
				rvals.push_back(255);
				gvals.push_back(0);
				bvals.push_back(0);
				camcoords.push_back(135);
			}
		}
		//travel to next gun  (64, 47, 70) to (33, 62, 77)
		start_x = 125;
		start_y = 109;
		start_z = 70;
		end_x = 250;
		end_y = 15;
		end_z = 77;
		total_duration = 0.035;
		start_time = time_coords.back();
		for (double i = 0; i < 4; i++)
		{

			time_coords.push_back(time_coords.back() + total_duration / 4);

			double interp = (double)(time_coords.back() - start_time) / total_duration;
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(0);
			bvals.push_back(0);
			camcoords.push_back(135);
		}
	}
	{
		//second shot
		double start_x = 250;
		double start_y = 15;
		double start_z = 70;
		double end_x = 125;
		double end_y = 109;
		double end_z = 70;
		double disruptor_percentage = 0.40;
		double total_duration = 0.249;
		double start_time = time_coords.back();
		//first shot
		for (double i = 0; i < total_duration; i += 1.0 / FPS)
		{

			double interp;

			time_coords.push_back(time_coords.back() + ((1.0 / FPS) / 2.0));

			// insert close point
			interp = max((1.0 + disruptor_percentage) / total_duration * (time_coords.back() - start_time) - disruptor_percentage, 0.0);
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(255);
			bvals.push_back(0);
			camcoords.push_back(135);

			// a little time has passed between points
			time_coords.push_back(time_coords.back() + ((1.0 / FPS) / 2.0));

			// insert far point
			interp = min((1.0 + disruptor_percentage) / total_duration * (time_coords.back() - start_time), 1.0);
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(255);
			bvals.push_back(0);
			camcoords.push_back(135);
		}
		//little explosion at (64, 47, 70) with -100 deg rotation
		start_x = 125;
		start_y = 109;
		start_z = 70;
		total_duration = 0.1295;
		start_time = time_coords.back();
		double rotation = 100.0 / 180.0 * PI;
		double scale = 4.0;
		for (double i = start_time; i < start_time + total_duration; i += 1.0 / FPS)
		{
			for (int j = 0; j < 8; j++)
			{
				time_coords.push_back(time_coords.back() + 1.0 / FPS / 8.0);

				xcoords.push_back(start_x + scale * (explosionxcoords[j % 8] * cos(rotation) + explosionycoords[j % 8] * sin(rotation)));
				ycoords.push_back(start_y + scale * (explosionycoords[j % 8] * cos(rotation) - explosionxcoords[j % 8] * sin(rotation)));
				zcoords.push_back(start_z);
				rvals.push_back(255);
				gvals.push_back(0);
				bvals.push_back(0);
				camcoords.push_back(135);
			}
		}
		//travel to next gun  (64, 47, 70) to (33, 62, 77)
		start_x = 125;
		start_y = 109;
		start_z = 70;
		end_x = 250;
		end_y = 15;
		end_z = 77;
		total_duration = 0.035;
		start_time = time_coords.back();
		for (double i = 0; i < 4; i++)
		{

			time_coords.push_back(time_coords.back() + total_duration / 4);

			double interp = (double)(time_coords.back() - start_time) / total_duration;
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(0);
			bvals.push_back(0);
			camcoords.push_back(135);
		}
	}
	{
		//second shot
		double start_x = 250;
		double start_y = 15;
		double start_z = 70;
		double end_x = 125;
		double end_y = 109;
		double end_z = 70;
		double disruptor_percentage = 0.40;
		double total_duration = 0.249;
		double start_time = time_coords.back();
		//first shot
		for (double i = 0; i < total_duration; i += 1.0 / FPS)
		{

			double interp;

			time_coords.push_back(time_coords.back() + ((1.0 / FPS) / 2.0));

			// insert close point
			interp = max((1.0 + disruptor_percentage) / total_duration * (time_coords.back() - start_time) - disruptor_percentage, 0.0);
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(255);
			bvals.push_back(0);
			camcoords.push_back(135);

			// a little time has passed between points
			time_coords.push_back(time_coords.back() + ((1.0 / FPS) / 2.0));

			// insert far point
			interp = min((1.0 + disruptor_percentage) / total_duration * (time_coords.back() - start_time), 1.0);
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(255);
			bvals.push_back(0);
			camcoords.push_back(135);
		}
		//little explosion at (64, 47, 70) with -100 deg rotation
		start_x = 125;
		start_y = 109;
		start_z = 70;
		total_duration = 0.1295;
		start_time = time_coords.back();
		double rotation = 100.0 / 180.0 * PI;
		double scale = 4.0;
		for (double i = start_time; i < start_time + total_duration; i += 1.0 / FPS)
		{
			for (int j = 0; j < 8; j++)
			{
				time_coords.push_back(time_coords.back() + 1.0 / FPS / 8.0);

				xcoords.push_back(start_x + scale * (explosionxcoords[j % 8] * cos(rotation) + explosionycoords[j % 8] * sin(rotation)));
				ycoords.push_back(start_y + scale * (explosionycoords[j % 8] * cos(rotation) - explosionxcoords[j % 8] * sin(rotation)));
				zcoords.push_back(start_z);
				rvals.push_back(255);
				gvals.push_back(0);
				bvals.push_back(0);
				camcoords.push_back(135);
			}
		}
	}

	//Shoot photon torpedo from (4,26,45) to (185,78,156), crossing through the point(113, 0, 196)
	//Should take 6 seconds @ 3.5ms/point = 1714 points
	{
		//travel to next gun  (64, 47, 70) to (4, 26, 45)
		double start_x = 125;
		double start_y = 109;
		double start_z = 70;
		double end_x = 0;
		double end_y = 16;
		double end_z = 45;
		double total_duration = 0.21;
		double start_time = time_coords.back();
		for (double i = 0; i < 4; i++)
		{

			time_coords.push_back(time_coords.back() + total_duration / 4);

			double interp = (double)(time_coords.back() - start_time) / total_duration;
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(0);
			bvals.push_back(0);
			camcoords.push_back(135);
		}

		// push in photon torpoedo
		// wait
		total_duration = 0.02;
		time_coords.push_back(time_coords.back() + total_duration);
		xcoords.push_back(end_x);
		ycoords.push_back(end_y);
		zcoords.push_back(end_z);
		rvals.push_back(0);
		gvals.push_back(0);
		bvals.push_back(255);
		camcoords.push_back(135);

		// second point
		total_duration = 0.84;
		time_coords.push_back(time_coords.back() + total_duration);
		xcoords.push_back(75);
		ycoords.push_back(75);
		zcoords.push_back(194);
		rvals.push_back(0);
		gvals.push_back(0);
		bvals.push_back(255);
		camcoords.push_back(135);

		// third point
		total_duration = 0.497;
		time_coords.push_back(time_coords.back() + total_duration);
		xcoords.push_back(150);
		ycoords.push_back(25);
		zcoords.push_back(156);
		rvals.push_back(0);
		gvals.push_back(0);
		bvals.push_back(255);
		camcoords.push_back(135);

		// fourth point
		total_duration = 0.2;
		time_coords.push_back(time_coords.back() + total_duration);
		xcoords.push_back(188);
		ycoords.push_back(109);
		zcoords.push_back(156);
		rvals.push_back(0);
		gvals.push_back(0);
		bvals.push_back(255);
		camcoords.push_back(135);
	
		//Big explosion!! at (185,78,156) with 90 deg rotation
		double rotation = -90.0 / 180.0 * PI;
		double scale;
		double pos_x = 188;
		double pos_y = 78;
		double pos_z = 156;
		total_duration = 1.484;

		for(double t = 0; t < total_duration;){
			for(int j = 0; j<10; j++){
				t += 1.0 / FPS / 10.0;
				double interp = (double)t / total_duration;
				scale = 4.0 * (1.0 - interp) + 9.0 * interp;

				time_coords.push_back(time_coords.back() + 1.0 / FPS / 10.0);
				xcoords.push_back(pos_x + scale * (explosion2xcoords[j % 10] * cos(rotation) + explosion2ycoords[j % 10] * sin(rotation)));
				ycoords.push_back(pos_y + scale * (explosion2ycoords[j % 10] * cos(rotation) - explosion2xcoords[j % 10] * sin(rotation)));
				zcoords.push_back(pos_z);
				rvals.push_back(255);
				gvals.push_back(0);
				bvals.push_back(0);
				camcoords.push_back(135);
			}
		}

		// std::cout << "Last point in x explosion: " << xcoords.back() << endl;
		// std::cout << "Last point in y explosion: " << ycoords.back() << endl;

		start_x = 225;
		start_y = 96;
		start_z = 156;
		end_x = 250;
		end_y = 15;
		end_z = 70;
		total_duration = 0.0315;
		start_time = time_coords.back();
		for (double i = 0; i < 4; i++)
		{

			time_coords.push_back(time_coords.back() + total_duration / 4);

			double interp = (double)(time_coords.back() - start_time) / total_duration;
			xcoords.push_back(start_x * (1 - interp) + end_x * interp);
			ycoords.push_back(start_y * (1 - interp) + end_y * interp);
			zcoords.push_back(start_z * (1 - interp) + end_z * interp);
			rvals.push_back(0);
			gvals.push_back(0);
			bvals.push_back(0);
			camcoords.push_back(135);
		}
		//one more at the end to make sure it has rested
		time_coords.push_back(time_coords.back() + total_duration / 4);
		xcoords.push_back(250);
		ycoords.push_back(15);
		zcoords.push_back(70);
		rvals.push_back(0);
		gvals.push_back(0);
		bvals.push_back(0);
		camcoords.push_back(135);
	}

	vector<double> x_interp;
	vector<double> y_interp;
	vector<double> z_interp;
	vector<double> r_interp;
	vector<double> g_interp;
	vector<double> b_interp;

	perform_interpolation(&time_coords, &xcoords, &x_interp);
	perform_interpolation(&time_coords, &ycoords, &y_interp);
	perform_interpolation(&time_coords, &zcoords, &z_interp);
	perform_interpolation(&time_coords, &rvals, &r_interp);
	perform_interpolation(&time_coords, &gvals, &g_interp);
	perform_interpolation(&time_coords, &bvals, &b_interp);

	clamp0to255(&x_interp);
	clamp0to255(&y_interp);
	clamp0to255(&z_interp);
	clamp0to255(&r_interp);
	clamp0to255(&g_interp);
	clamp0to255(&b_interp);

	//output x coords
	/*std::cout << "X_interp" << endl;
	for(size_t i = 0; i < x_interp.size(); i++){
		std::cout << +static_cast<uint8_t>(x_interp.at(i)) << "\t";
	}*/

	//output x coords
	std::cout << endl
			  << endl
			  << "X" << endl;
	for (size_t i = 0; i < x_interp.size(); i++)
	{
		std::cout << "," << +static_cast<uint8_t>(x_interp.at(i) + .5);
	}

	//output y coords
	std::cout << endl
			  << endl
			  << "Y" << endl;
	for (size_t i = 0; i < y_interp.size(); i++)
	{
		std::cout << "," << +static_cast<uint8_t>(255.0 - y_interp.at(i) + .5);
	}

	//output z coords
	std::cout << endl
			  << endl
			  << "Z" << endl;
	for (size_t i = 0; i < z_interp.size(); i++)
	{
		std::cout << "," << +static_cast<uint8_t>(z_interp.at(i) + .5);
	}

	//output r coords
	std::cout << endl
			  << endl
			  << "R" << endl;
	for (size_t i = 0; i < r_interp.size(); i++)
	{
		std::cout << "," << +static_cast<uint8_t>(r_interp.at(i) + .5);
	}

	//output g coords
	std::cout << endl
			  << endl
			  << "G" << endl;
	for (size_t i = 0; i < g_interp.size(); i++)
	{
		std::cout << "," << +static_cast<uint8_t>(g_interp.at(i) + .5);
	}

	//output b coords
	std::cout << endl
			  << endl
			  << "B" << endl;
	for (size_t i = 0; i < b_interp.size(); i++)
	{
		std::cout << "," << +static_cast<uint8_t>(b_interp.at(i) + .5);
	}

	//output time coords
	std::cout << endl
			  << endl
			  << "TimeCoords" << endl;
	for (size_t i = 0; i < time_coords.size(); i++)
	{
		std::cout << time_coords.at(i) << "\t";
	}

	/*camcoords = {};
	for(int i = 0; i < 1000; i++){
		double t = (double)i / 1000;
		camcoords.push_back(0 * (2 * t * t * t - 3 * t * t + 1) + 0 * (t * t * t - 2 * t * t + t) + 173 * (-2 * t * t * t + 3 * t * t) + 0 * (t * t * t - t * t));
	}
	for(int i = 0; i < 606; i++){
		double t = (double)i / 606;
		camcoords.push_back(173 * (2 * t * t * t - 3 * t * t + 1) + 0 * (t * t * t - 2 * t * t + t) + 55 * (-2 * t * t * t + 3 * t * t) + -200 * (t * t * t - t * t));
	}
	for(int i = 0; i < 434; i++){
		double t = (double)i / 434;
		camcoords.push_back(55 * (2 * t * t * t - 3 * t * t + 1) + 1000 * (t * t * t - 2 * t * t + t) + 0 * (-2 * t * t * t + 3 * t * t) + 0 * (t * t * t - t * t));
	}*/

	//output the size of the arrays
	std::cout << endl;
	std::cout << "Size of arrays: " << b_interp.size();

	//output cam coords
	std::cout << endl
			  << endl
			  << "CAM" << endl;
	for (size_t i = 0; i < camcoords.size(); i++)
	{
		std::cout << "," << camcoords.at(i);
	}
	//output r coords
	std::cout << endl
			  << endl
			  << "R2" << endl;
	for (size_t i = 0; i < rvals.size(); i++)
	{
		std::cout << "," << rvals.at(i);
	}
	//output g coords
	std::cout << endl
			  << endl
			  << "G2" << endl;
	for (size_t i = 0; i < gvals.size(); i++)
	{
		std::cout << "," << gvals.at(i);
	}
	std::cout << endl
			  << endl
			  << camcoords.size();
	system("pause");
}
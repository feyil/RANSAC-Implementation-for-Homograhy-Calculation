// 230201057

// ------------------------------
// Written by Mustafa Ozuysal
// Contact <mustafaozuysal@iyte.edu.tr> for comments and bug reports
// ------------------------------
// Copyright (c) 2018, Mustafa Ozuysal
// All rights reserved.
// ------------------------------
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the copyright holders nor the
//       names of his/its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// ------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ------------------------------
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "image.hpp"
#include "homography.hpp"

using ceng391::Image;
using ceng391::short_to_image;
using ceng391::Keypoint;
using ceng391::Descriptor;
using ceng391::Match;
using ceng391::fit_homography4;
using std::vector;
using std::cout;
using std::endl;
using std::printf;

int main(int argc, char** argv)
{
     Image wall1(4, 4, 1);
     Image wall2(4, 4, 1);
     wall1.read_pnm("../wall1.ppm");
     wall2.read_pnm("../wall2.ppm");

     // double theta = 0 * 3.1415926 / 180;
     // Image wall2(wall3.w() + 100, wall3.h() + 100, 1);
     // wall3.rotate_centered(&wall2, theta);

     wall1.to_grayscale();
     wall2.to_grayscale();

     float threshold = 1e6;
     float k = 0.06f;
     float sigma = 2.5f;
     vector<Keypoint> keys1 = wall1.harris_corners(threshold, k, sigma);
     vector<Keypoint> keys2 = wall2.harris_corners(threshold, k, sigma);
     cout << "Detected " << keys1.size()
          << " keypoints on wall1 and " << keys2.size()
          << " keypoints on wall2" << endl;

     vector<Descriptor> desc1 = wall1.compute_brief(keys1);
     vector<Descriptor> desc2 = wall2.compute_brief(keys2);
     cout << "Extracted " << desc1.size()
          << " descriptors on wall1 and " << desc2.size()
          << " descriptors on wall2" << endl;

     vector<Match> wall_matches = Image::match_brief(desc1, desc2);

     // Calculating homography
     double h[9] = {};
     int n_inliers_best = Image::find_homography(wall_matches, keys1, keys2, h);

     for (int i = 0; i < 9; ++i) {
               h[i] /= h[8];
     }
     cout<<"n_inliers_best = "<<n_inliers_best<<endl;
     printf("---------- H ----------\n");
     printf("%6.2f %6.2f %6.2f\n", h[0], h[3], h[6]);
     printf("%6.2f %6.2f %6.2f\n", h[1], h[4], h[7]);
     printf("%6.2f %6.2f %6.2f\n", h[2], h[5], h[8]);

     // Transform with the best homography and select inliers
     vector<Match> ransac_matches;
     for(int i = 0; i < wall_matches.size(); i++) {

          Keypoint keypoint1 = keys1[wall_matches[i].key_id0];
          Keypoint keypoint2 = keys2[wall_matches[i].key_id1];

          double w_h = h[2] * keypoint1.x + h[5] * keypoint1.y + h[8];
          double x_h = (h[0] * keypoint1.x + h[3] * keypoint1.y + h[6]) / w_h;
          double y_h = (h[1] * keypoint1.x + h[4] * keypoint1.y + h[7]) / w_h;

          if(abs(keypoint2.x - x_h) < 3 && abs(keypoint2.y - y_h) < 3) {
               
               ransac_matches.push_back(wall_matches[i]);

               // cout<<"Match Number: "<<i<<endl;
               // cout<<"\tkey_id0xH:"<<"("<<x_h<<","<<y_h<<")"<<endl;
               // cout<<"\tkey_id1:"<<"("<<keypoint2.x<<","<<keypoint2.y<<")"<<endl<<endl;
          }
     }

     // Looking how well ransac is performing its job
     Image *match_image = make_match_image(wall_matches, wall1, keys1,
                                              wall2, keys2, 2000);
     match_image->write_pnm("/tmp/wall_matches");
     

     match_image = make_match_image(ransac_matches, wall1, keys1,
                                             wall2, keys2, 2000);
     match_image->write_pnm("/tmp/wall_matches_ransac");

     cout<<"Image written to /tmp/wall_matches"<<endl;
     cout<<"Image written to /tmp/wall_matches_ransac"<<endl;

     delete match_image;

     return EXIT_SUCCESS;
}


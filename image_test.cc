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
     //    Image* gray = Image::new_gray(128, 128);
     //    cout << "(" << gray->w() << "x" << gray->h() << ") channels: "
     //         << gray->n_ch() << " step: " << gray->step() << endl;
     //    gray->set_zero();
     //    gray->set_rect(32, 32, 64, 64, 255);
     //    gray->write_pnm("/tmp/test_image_gray");
     //    delete gray;

     //    Image* rgb = Image::new_rgb(128, 128);
     //    cout << "(" << rgb->w() << "x" << rgb->h() << ") channels: "
     //         << rgb->n_ch() << " step: " << rgb->step() << endl;
     //    rgb->set_zero();
     //    rgb->set_rect(32, 32, 64, 64, 255, 0, 255);
     //    rgb->write_pnm("/tmp/test_image_rgb");
     //    delete rgb;

     //    Image img(4, 4, 1);
     //    img.read_pnm("/tmp/test_image_gray.pgm");
     //    img.to_rgb();
     //    img.write_pnm("/tmp/test_image_gray2rgb");

     //    img.read_pnm("/tmp/test_image_rgb.ppm");
     //    img.to_grayscale();
     //    img.write_pnm("/tmp/test_image_rgb2gray");

     //    img.read_pnm("../small_city.pgm");
     //    Image rotated(img.w()*2, img.h()*2, 1);
     //    double theta = 45.0 * 3.1415926 / 180;
     //    img.rotate_centered(&rotated, theta);
     //    rotated.write_pnm("/tmp/small_city_crotated_45");

     //    float threshold = 1000.0f;
     //    float k = 0.06f;
     //    float sigma = 2.5f;
     //    vector<Keypoint> keys = img.harris_corners(threshold, k, sigma);
     //    cout << "Detected " << keys.size()
     //         << " keypoints on small_city.pgm" << endl;
     //    Image *key_image = make_keypoint_image(&img, &keys);
     //    key_image->write_pnm("/tmp/keys");

     //    short *dx = img.deriv_x();
     //    short *dy = img.deriv_y();
     //    cout << "Derivatives computed" << endl;

     //    Image *dx_img = short_to_image(dx, img.w(), img.h());
     //    Image *dy_img = short_to_image(dy, img.w(), img.h());
     //    dx_img->write_pnm("/tmp/dx");
     //    dy_img->write_pnm("/tmp/dy");

     //    float sigma_x = 5.5f;
     //    float sigma_y = 5.5f;
     //    img.smooth(sigma_x, sigma_y);
     //    img.write_pnm("/tmp/smoothed_xy");

     Image wall1(4, 4, 1);
     Image wall2(4, 4, 1);
     wall1.read_pnm("../wall1.ppm");
     wall2.read_pnm("../wall2.ppm");
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
     Image *match_image = make_match_image(wall_matches, wall1, keys1,
                                              wall2, keys2, 20);
     match_image->write_pnm("/tmp/wall_matches");

     double matches[16] = {
               0.0, 0.0, 3.0, 2.0,
               1.0, 0.0, 4.0, 2.0,
               1.0, 1.0, 4.0, 3.0,
               0.0, 1.0, 3.0, 3.0
     };
     double h[9];
     fit_homography4(&matches[0], &h[0]);
     for (int i = 0; i < 9; ++i) {
               h[i] /= h[8];
     }
     printf("---------- H ----------\n");
     printf("%6.2f %6.2f %6.2f\n", h[0], h[3], h[6]);
     printf("%6.2f %6.2f %6.2f\n", h[1], h[4], h[7]);
     printf("%6.2f %6.2f %6.2f\n", h[2], h[5], h[8]);

     //    delete key_image;
     //    delete [] dx;
     //    delete [] dy;
     //    delete dx_img;
     //    delete dy_img;

     delete match_image;

     return EXIT_SUCCESS;
}

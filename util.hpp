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
#ifndef UTIL_HPP
#define UTIL_HPP

#include <vector>

#include "types.hpp"
#include "keypoint.hpp"
#include "descriptor.hpp"

namespace ceng391 {

class Image;

Image *short_to_image(const short *ptr, int width, int height);

float *gaussian_kernel(float sigma, int *k);

template <typename T>
void copy_to_buffer(float *buffer, const T *src, int n,
                    int border_size, int stride);

template <typename T>
void copy_from_buffer(T *dest, const float *buffer, int n, int stride);

void convolve_buffer(float *buffer, int n, const float *kernel, int k);


short *vec_mul(int n, const short *v0, const short *v1);
void smooth_short_buffer(int w, int h, short *I, float sigma);
float *harris_corner_score(int w, int h, const short *Ix2, const short *Iy2,
                           const short *IxIy, float k);
Image *make_keypoint_image(Image *img, std::vector<Keypoint> *keys);
Image *make_match_image(const std::vector<Match> &matches,
                        const Image &img0, const std::vector<Keypoint> &keys0,
                        const Image &img1, const std::vector<Keypoint> &keys1,
                        int distance_threshold);
// ---------------------------- Template Definitions ------------------------------
template <typename T>
void copy_to_buffer(float *buffer, const T *src, int n, int border_size,
                    int stride)
{
        for (int i = 0; i < n; ++i)
                buffer[border_size + i] = src[i * stride];

        for (int i = 0; i < border_size; ++i) {
                buffer[i] = buffer[border_size];
                buffer[i + n + border_size] = buffer[border_size + n - 1];
        }
}

template <typename T>
void copy_from_buffer(T *dest, const float *buffer, int n, int stride)
{
        for (int i = 0; i < n; ++i)
                dest[i * stride] = buffer[i];
}

}

#endif

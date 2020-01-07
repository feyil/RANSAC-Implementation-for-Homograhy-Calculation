// ------------------------------
// Written by Mustafa Ozuysal
// Contact <mustafaozuysal@iyte.edu.tr> for comments and bug reports
// ------------------------------
// Copyright (c) 2019, Mustafa Ozuysal
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
#include "util.hpp"

#include <cstdio>
#include <cmath>
#include <memory>

#include "image.hpp"

using std::printf;
using std::max;
using std::ceil;
using std::exp;
using std::unique_ptr;

namespace ceng391 {

Image *short_to_image(const short *ptr, int width, int height)
{
        Image *img = Image::new_gray(width, height);
        for (int y = 0; y < height; ++y) {
                const short *row = ptr + y * width;
                uchar *irow = img->data(y);
                for (int x = 0; x < width; ++x) {
                        irow[x] = (row[x] + 255) / 2;
                }
        }

        return img;
}

float *gaussian_kernel(float sigma, int *k)
{
        int l = ceil(2.0f * sigma);
        *k = 2 * l + 1;
        float *kernel = new float[*k];
        float sum = 0.0f;
        for (int i = 0; i < *k; ++i) {
                int x = i - l;
                kernel[i] = exp(-0.5f * x * x / sigma / sigma);
                sum += kernel[i];
        }
        for (int i = 0; i < *k; ++i)
                kernel[i] /= sum;

        return kernel;
}

void convolve_buffer(float *buffer, int n, const float *kernel, int k)
{
        for (int i = 0; i < n; ++i) {
                float sum = 0.0f;
                for (int j = 0; j < k; ++j) {
                        sum += kernel[j] * buffer[i + j];
                }
                buffer[i] = sum;
        }
}

short *vec_mul(int n, const short *v0, const short *v1)
{
        short *p = new short[n];
        for (int i = 0; i < n; ++i) {
                p[i] = v0[i] * v1[i];
        }
        return p;
}

void smooth_short_buffer(int w, int h, short *I, float sigma)
{
        int k = 0;
        unique_ptr<float []> kernel(gaussian_kernel(sigma, &k));

        int l = k / 2;
        int max_wh = max(w, h);
        unique_ptr<float []>  buffer(new float[max_wh + 2 * l]);

        for (int y = 0; y < h - 1; ++y) {
                copy_to_buffer<short>(buffer.get(), I + y * w, w, l, 1);
                convolve_buffer(buffer.get(), w, kernel.get(), k);
                copy_from_buffer<short>(I + y * w, buffer.get(), w, 1);
        }

        for (int x = 0; x < w - 1; ++x) {
                copy_to_buffer<short>(buffer.get(), I + x, h, l, w);
                convolve_buffer(buffer.get(), h, kernel.get(), k);
                copy_from_buffer<short>(I + x, buffer.get(), h, w);
        }
}

float *harris_corner_score(int w, int h, const short *Ix2, const short *Iy2,
                           const short *IxIy, float k)
{
        float *score = new float[w * h];
        for (int y = 0; y < h; ++y) {
                const short *A = Ix2 + y * w;
                const short *B = IxIy + y * w;
                const short *C = Iy2 + y * w;
                float *R = score + y * w;
                for (int x = 0; x < w; ++x) {
                        float det = A[x] * C[x] - B[x] * B[x];
                        float tr = A[x] + C[x];
                        R[x] = det - k * tr * tr;
                }
        }

        return score;
}

static inline void paint_pixel(Image *img, int x, int y, uchar *color)
{
        img->data(y)[3*x] = color[0];
        img->data(y)[3*x + 1] = color[1];
        img->data(y)[3*x + 2] = color[2];
}

Image *make_keypoint_image(Image *img, std::vector<Keypoint> *keys)
{
        Image *rgb = Image::new_copy(img);
        rgb->to_rgb();

        uchar color[3] = { 255, 0, 0 };
        for (size_t i = 0; i < keys->size(); ++i) {
                int x = (*keys)[i].x;
                int y = (*keys)[i].y;

                if (x < 2 || x >= rgb->w() - 2
                    || y < 2 || y >= rgb->h()) {
                        continue;
                }

                paint_pixel(rgb, x - 2, y, &color[0]);
                paint_pixel(rgb, x + 2, y, &color[0]);
                paint_pixel(rgb, x - 1, y - 1, &color[0]);
                paint_pixel(rgb, x + 1, y - 1, &color[0]);
                paint_pixel(rgb, x, y - 2, &color[0]);
                paint_pixel(rgb, x - 1, y + 1, &color[0]);
                paint_pixel(rgb, x + 1, y + 1, &color[0]);
                paint_pixel(rgb, x, y + 2, &color[0]);
        }

        return rgb;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

static inline void paint_line(Image *img, int x0, int y0,
                              int x1, int y1, uchar *color)
{
        float deltax = x1 - x0;
        float deltay = y1 - y0;
        float deltaerr = fabs(deltay / deltax);
        float error = 0.0f;
        int y = y0;
        for (int x = x0; x <= x1; ++x) {
                paint_pixel(img, x, y, color);
                error += deltaerr;
                if (error >= 0.5) {
                        y += sgn(deltay) * 1.0f;
                        error -= 1.0f;
                }
        }
}

Image *make_match_image(const std::vector<Match> &matches,
                        const Image &img0, const std::vector<Keypoint> &keys0,
                        const Image &img1, const std::vector<Keypoint> &keys1,
                        int distance_threshold)
{
        if (img0.n_ch() != 1 || img1.n_ch() != 1) {
                fprintf(stderr, "Match image needs grayscale input!\n");
                exit(9);
        }

        int w = img0.w() + img1.w();
        int h = max(img0.h(), img1.h());
        Image *rgb = Image::new_rgb(w, h);

        rgb->set_zero();
        for (int y = 0; y < img0.h(); ++y) {
                uchar *crow = rgb->data(y);
                const uchar *row0 = img0.data(y);
                for (int x = 0; x < img0.w(); ++x) {
                        crow[3 * x + 0] = row0[x];
                        crow[3 * x + 1] = row0[x];
                        crow[3 * x + 2] = row0[x];
                }
        }

        for (int y = 0; y < img1.h(); ++y) {
                uchar *crow = rgb->data(y);
                const uchar *row1 = img1.data(y);
                for (int x = 0; x < img1.w(); ++x) {
                        crow[3 * (x + img0.w()) + 0] = row1[x];
                        crow[3 * (x + img0.w()) + 1] = row1[x];
                        crow[3 * (x + img0.w()) + 2] = row1[x];
                }
        }

        uchar color[3] = { 255, 0, 0 };
        for (size_t i = 0; i < matches.size(); ++i) {
                if (matches[i].distance > distance_threshold)
                        continue;

                int k0 = matches[i].key_id0;
                int k1 = matches[i].key_id1;
                int x0 = keys0[k0].x;
                int y0 = keys0[k0].y;
                int x1 = keys1[k1].x + img0.w();
                int y1 = keys1[k1].y;

                paint_line(rgb, x0, y0, x1, y1, &color[0]);
        }

        return rgb;
}


}

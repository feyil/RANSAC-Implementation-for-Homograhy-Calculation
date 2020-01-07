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
#include "image.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <ostream>
#include <memory>
#include <limits>

#include "brief_data.hpp"
#include "util.hpp"
#include "homography.hpp"

using std::cerr;
using std::clog;
using std::cos;
using std::min;
using std::endl;
using std::exp;
using std::ifstream;
using std::ios;
using std::memset;
using std::ofstream;
using std::sin;
using std::string;
using std::unique_ptr;
using std::numeric_limits;

namespace ceng391 {

Image::Image(int width, int height, int n_channels, int step)
{
        m_width = width;
        m_height = height;
        m_n_channels = n_channels;

        m_step = m_width*m_n_channels;
        if (m_step < step)
                m_step = step;
        m_data = new uchar[m_step*height];
}

Image::~Image()
{
        delete [] m_data;
}

Image* Image::new_gray(int width, int height)
{
        return new Image(width, height, 1);
}

Image* Image::new_rgb(int width, int height)
{
        return new Image(width, height, 3);
}

Image *Image::new_copy(Image *img)
{
        Image *cpy = new Image(img->w(), img->h(), img->n_ch());
        for (int y = 0; y < img->h(); ++y)
                memcpy(cpy->data(y), img->data(y), img->w() * img->n_ch());
        return cpy;
}

void Image::set_rect(int x, int y, int width, int height, uchar red, uchar green, uchar blue)
{
        if (x < 0) {
                width += x;
                x = 0;
        }

        if (y < 0) {
                height += y;
                y = 0;
        }

        if (m_n_channels == 1) {
                int value = 0.3*red + 0.59*green + 0.11*blue;
                if (value > 255)
                        value = 255;
                for (int j = y; j < y+height; ++j) {
                        if (j >= m_height)
                                break;
                        uchar* row_data = data(j);
                        for (int i = x; i < x+width; ++i) {
                                if (i >= m_width)
                                        break;
                                row_data[i] = value;
                        }
                }
        } else if (m_n_channels == 3) {
                for (int j = y; j < y+height; ++j) {
                        if (j >= m_height)
                                break;
                        uchar* row_data = data(j);
                        for (int i = x; i < x+width; ++i) {
                                if (i >= m_width)
                                        break;
                                row_data[i*3]     = red;
                                row_data[i*3 + 1] = green;
                                row_data[i*3 + 2] = blue;
                        }
                }
        }
}

void Image::set_rect(int x, int y, int width, int height, uchar value)
{
        if (x < 0) {
                width += x;
                x = 0;
        }

        if (y < 0) {
                height += y;
                y = 0;
        }

        for (int j = y; j < y+height; ++j) {
                if (j >= m_height)
                        break;
                uchar* row_data = data(j);
                for (int i = x; i < x+width; ++i) {
                        if (i >= m_width)
                                break;
                        for (int c = 0; c < m_n_channels; ++c)
                                row_data[i*m_n_channels + c] = value;
                }
        }
}

void Image::to_grayscale()
{
        if (m_n_channels == 1) {
                return;
        } else if (m_n_channels == 3) {
                int new_step = m_width;
                uchar *new_data = new uchar[new_step * m_height];
                for (int y = 0; y < m_height; ++y) {
                        uchar *row_old = m_data + m_step * y;
                        uchar *row_new = new_data + new_step * y;
                        for (int x = 0; x < m_width; ++x) {
                                uchar red = row_old[3*x];
                                uchar green = row_old[3*x + 1];
                                uchar blue = row_old[3*x + 2];
                                int value = 0.3*red + 0.59*green + 0.11*blue;
                                if (value > 255)
                                        value = 255;
                                row_new[x] = value;
                        }
                }

                delete [] m_data;
                m_data = new_data;
                m_step = new_step;
                m_n_channels = 1;
        }
}

void Image::to_rgb()
{
        if (m_n_channels == 3) {
                return;
        } else if (m_n_channels == 1) {
                int new_step = m_width * 3;
                uchar *new_data = new uchar[new_step * m_height];
                for (int y = 0; y < m_height; ++y) {
                        uchar *row_old = m_data + m_step * y;
                        uchar *row_new = new_data + new_step * y;
                        for (int x = 0; x < m_width; ++x) {
                                uchar value = row_old[x];
                                row_new[3*x]     = value;
                                row_new[3*x + 1] = value;
                                row_new[3*x + 2] = value;
                        }
                }

                delete [] m_data;
                m_data = new_data;
                m_step = new_step;
                m_n_channels = 3;
        }
}

bool Image::write_pnm(const std::string& filename) const
{
        ofstream fout;

        string magic_head = "P5";
        string extended_name = filename + ".pgm";
        if (m_n_channels == 3) {
                magic_head = "P6";
                extended_name = filename + ".ppm";
        }

        fout.open(extended_name.c_str(), ios::out | ios::binary);
        if (!fout.good()) {
                cerr << "Error opening file " << extended_name << " for output!" << endl;
                return false;
        }

        fout << magic_head << "\n";
        fout << m_width << " " << m_height << " 255\n";
        for (int y = 0; y < m_height; ++y) {
                const uchar *row_data = data(y);
                fout.write(reinterpret_cast<const char*>(row_data), m_width*m_n_channels*sizeof(uchar));
        }
        fout.close();

        return true;
}

bool Image::read_pnm(const std::string& filename)
{
        ifstream fin(filename.c_str(), ios::in | ios::binary);
        if (!fin.good()) {
                cerr << "Error opening PNM file " << filename << endl;
                return false;
        }

        int width;
        int height;
        int max_val;
        int n_channels = 1;
        string head = "00";
        head[0] = fin.get();
        head[1] = fin.get();
        if (head == "P5") {
                clog << "Loading PGM Binary" << endl;
                n_channels = 1;
        } else if (head == "P6") {
                clog << "Loading PPM Binary" << endl;
                n_channels = 3;
        } else {
                cerr << "File " << filename << " is not a Binary PGM or PPM!" << endl;
                return false;
        }

        fin >> width;
        fin >> height;
        fin >> max_val;
        if (fin.peek() == '\n')
                fin.get();

        int step = width * n_channels;
        uchar *new_data = new uchar[step*height];
        for (int y = 0; y < height; ++y) {
                fin.read(reinterpret_cast<char*>(new_data + y*step), step*sizeof(uchar));
        }
        fin.close();

        delete [] m_data;
        m_data = new_data;
        m_width = width;
        m_height = height;
        m_step = step;
        m_n_channels = n_channels;

        return true;
}

short *Image::deriv_x() const
{
        if (m_n_channels == 3) {
                cerr << "Image derivatives only implemented for grayscale images!" << endl;
                return nullptr;
        }

        short *dx = new short[m_width * m_height];
        for (int y = 0; y < m_height; ++y) {
                const uchar *row = this->data(y);
                short *drow = dx + y * m_width;
                drow[0] = 0;
                for (int x = 1; x < m_width - 1; ++x) {
                        drow[x] = row[x + 1] - row[x - 1];
                }
                drow[m_width - 1] = 0;
        }

        return dx;
}

short *Image::deriv_y() const
{
        if (m_n_channels == 3) {
                cerr << "Image derivatives only implemented for grayscale images!" << endl;
                return nullptr;
        }

        short *dy = new short[m_width * m_height];

        memset(dy, 0, m_width * sizeof(*dy));
        for (int y = 1; y < m_height - 1; ++y) {
                const uchar *rowm = this->data(y - 1);
                const uchar *rowp = this->data(y + 1);
                short *drow = dy + y * m_width;
                for (int x = 0; x < m_width; ++x) {
                        drow[x] = rowp[x] - rowm[x];
                }
        }
        memset(dy + (m_height - 1) * m_width, 0, m_width * sizeof(*dy));

        return dy;
}

void Image::rotate(Image *rotated, double theta, double tx, double ty) const
{
        if (m_n_channels != 1) {
                cerr << "Rotate only works on grayscale images!" << endl;
                return;
        }
        rotated->to_grayscale();

        double ct = cos(theta);
        double st = sin(theta);
        double tx_inv = -ct * tx + st * ty;
        double ty_inv = -st * tx - ct * ty;

        int wp = rotated->w();
        int hp = rotated->h();

        for (int yp = 0; yp < hp; ++yp) {
                uchar *row_p = rotated->data(yp);
                for (int xp = 0; xp < wp; ++xp) {
                        double x = ct * xp - st * yp + tx_inv;
                        double y = st * xp + ct * yp + ty_inv;

                        int x0 = static_cast<int>(x);
                        int y0 = static_cast<int>(y);

                        int value = 0;
                        if (x0 < 0 || y0 < 0 || x0 >= m_width || y0 >= m_height) {
                                value = 0;
                        } else {
                                const uchar *row = this->data(y0);
                                value = row[x0];
                        }

                        row_p[xp] = value;
                }
        }
}

void Image::rotate_centered(Image *rotated, double theta) const
{
        double ct = cos(theta);
        double st = sin(theta);
        double hw = m_width / 2.0;
        double hh = m_height / 2.0;
        double hwp = rotated->w() / 2.0;
        double hhp = rotated->h() / 2.0;

        double tx_cap = -ct * hw - st * hh + hwp;
        double ty_cap =  st * hw - ct * hh + hhp;
        this->rotate(rotated, theta, tx_cap, ty_cap);
}

void Image::smooth_x(float sigma)
{
        if (m_n_channels != 1) {
                cerr << "Smooth-x only works on grayscale images!" << endl;
                return;
        }

        int k = 0;
        unique_ptr<float []> kernel(gaussian_kernel(sigma, &k));

        int l = k / 2;
        unique_ptr<float []>  buffer(new float[m_width + 2 * l]);

        for (int y = 0; y < m_height - 1; ++y) {
                copy_to_buffer(buffer.get(), this->data(y), m_width, l, 1);
                convolve_buffer(buffer.get(), m_width, kernel.get(), k);
                copy_from_buffer(this->data(y), buffer.get(), m_width, 1);
        }
}

void Image::smooth_y(float sigma)
{
        if (m_n_channels != 1) {
                cerr << "Smooth-x only works on grayscale images!" << endl;
                return;
        }

        int k = 0;
        unique_ptr<float []> kernel(gaussian_kernel(sigma, &k));

        int l = k / 2;
        unique_ptr<float []>  buffer(new float[m_height + 2 * l]);

        for (int x = 0; x < m_width - 1; ++x) {
                copy_to_buffer(buffer.get(), m_data + x, m_height, l, m_step);
                convolve_buffer(buffer.get(), m_height, kernel.get(), k);
                copy_from_buffer(m_data + x, buffer.get(), m_height, m_step);
        }
}

void Image::smooth(float sigma_x, float sigma_y)
{
        smooth_x(sigma_x);
        smooth_y(sigma_y);
}

std::vector<Keypoint> Image::harris_corners(float threshold, float k,
                                            float sigma)
{
        short *Ix = this->deriv_x();
        short *Iy = this->deriv_y();

        int n = m_width * m_height;
        short *Ix2  = vec_mul(n, Ix, Ix);
        short *Iy2  = vec_mul(n, Iy, Iy);
        short *IxIy = vec_mul(n, Ix, Iy);

        smooth_short_buffer(m_width, m_height, Ix2, sigma);
        smooth_short_buffer(m_width, m_height, Iy2, sigma);
        smooth_short_buffer(m_width, m_height, IxIy, sigma);

        const float *score = harris_corner_score(m_width, m_height,
                                                 Ix2, Iy2, IxIy, k);

        std::vector<Keypoint> keys;
        int border = min(1, static_cast<int>(2.0f * sigma));
        for (int y = border; y < m_height - border; ++y) {
                const float *Rm = score + (y - 1) * m_width;
                const float *R  = score + y * m_width;
                const float *Rp = score + (y + 1) * m_width;
                for (int x = border; x < m_width - border; ++x) {
                        if (R[x] > threshold && R[x] > R[x - 1]
                            && R[x] > R[x + 1] && R[x] > Rm[x]
                            && R[x] > Rp[x]) {
                                Keypoint key;
                                key.x = x;
                                key.y = y;
                                key.score = R[x];
                                keys.push_back(key);
                        }
                }
        }

        delete [] score;
        delete [] IxIy;
        delete [] Iy2;
        delete [] Ix2;
        delete [] Iy;
        delete [] Ix;

        return keys;
}

std::vector<Descriptor> Image::compute_brief(const std::vector<Keypoint>& keys)
{
        Image *cpy = Image::new_copy(this);
        cpy->to_grayscale();
        cpy->smooth(2.5f, 2.5f);

        int n = static_cast<int>(keys.size());
        vector<Descriptor> descriptors;
        descriptors.reserve(n); // preallocates space for n descriptors, does not change size

        for (int i = 0; i < n; ++i) {
                int x = keys[i].x;
                int y = keys[i].y;
                if (x < 8 || x >= (cpy->w() - 8)
                    || y < 8 || (y >= cpy->h() - 8)) {
                        continue; // skip keypoint too close to border
                }

                Descriptor d;
                d.key_id = i;
                for (int j = 0; j < 32; ++j) {
                        uchar value = 0;
                        for (int k = 0; k < 8; ++k) {
                                int bit_index = j * 8 + k;
                                int x0 = x + BRIEF_OFFSETS[bit_index][0][0];
                                int y0 = y + BRIEF_OFFSETS[bit_index][0][1];
                                int x1 = x + BRIEF_OFFSETS[bit_index][1][0];
                                int y1 = y + BRIEF_OFFSETS[bit_index][1][1];
                                uchar I0 = cpy->data(y0)[x0];
                                uchar I1 = cpy->data(y1)[x1];
                                if (I0 > I1)
                                        value |= (1 << k);
                        }
                        d.desc[j] = value;
                }
                descriptors.push_back(d);
        }

        delete cpy;

        return descriptors;
}


// From Hacker's Delight, 2nd Ed., Section 5.1.  You could improve this quite a
// bit from by using the information in the book.
static inline int popcnt8(uchar x)
{
        x = (x & 0x55) + ((x >> 1) & 0x55);
        x = (x & 0x33) + ((x >> 2) & 0x33);
        x = (x & 0x0F) + (x >> 4);

        return x;
}

static inline int hamming_distance(int n, const uchar *v0, const uchar *v1)
{
        int d = 0;
        for (int i = 0; i < n; ++i) {
                uchar x = v0[i] ^ v1[i];
                d += popcnt8(x);
        }

        return d;
}

std::vector<Match> Image::match_brief(const std::vector<Descriptor>& d0,
                                      const std::vector<Descriptor>& d1)
{
        int n0 = static_cast<int>(d0.size());
        int n1 = static_cast<int>(d1.size());

        vector<Match> matches;
        matches.reserve(n0); // preallocates space for n0 matches
        for (int i = 0; i < n0; ++i) {
                int best_id = -1;
                int best_distance = numeric_limits<int>::max();

                for (int j = 0; j < n1; ++j) {
                        int d = hamming_distance(32, &d0[i].desc[0],
                                                 &d1[j].desc[0]);
                        if (d < best_distance) {
                                best_distance = d;
                                best_id = j;
                        }
                }

                Match m;
                m.key_id0 = d0[i].key_id;
                m.key_id1 = d1[best_id].key_id;
                m.distance = best_distance;
                matches.push_back(m);
        }

        return matches;
}

int Image::find_homography(const std::vector<Match> &matches, 
                    const std::vector<Keypoint> &keys1,
                    const std::vector<Keypoint> &keys2,
                    double* h)
{
        int n_inliers_best = 0;
        double h_best[9] = { }; //Initialized to 0
        
        const int FIXED_ITERATION_SIZE = 1;
        
        for(int i = 0; i < FIXED_ITERATION_SIZE; i++) {
                
                // Finding random four candidate match for homography calculation
                double candidate_matches[16];
                for(int j = 0; j < 4; j++) {
                        int candidate = generate_random(0, matches.size() - 1);

                        Keypoint keypoint1 = keys1[matches[candidate].key_id0];
                        Keypoint keypoint2 = keys2[matches[candidate].key_id1];

                        candidate_matches[4 * j] = keypoint1.x;
                        candidate_matches[4 * j + 1] = keypoint1.y;

                        candidate_matches[4 * j + 2] = keypoint2.x;
                        candidate_matches[4* j + 3] = keypoint2.y;
                }

                
                
                

                // use the fit_homography4 to compute a homography that fits the four selected matches

                // transform with homography

                // is_inlier

                // update the inlier

        }

        return n_inliers_best;
}

}

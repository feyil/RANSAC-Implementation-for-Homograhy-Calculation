# RANSAC Algorihtm Implementation

### Introduction

* This repo includes my solution of the given homework(5/5) in the scope of the Introduction to Image Understanding(CENG391) course which is given as a technical elective in 2019-2020 Fall semester by Computer Engineering Department at Izmir Institute of Technology.

* (*)README.md file uses some parts of the official Homework Doc to better express the purpose of the Homework.

* All solutions implemented on top of the base code.


### Problem*

* Write a new member function Image::find_homography that takes a const reference to a std::vector of matches, two const references to std::vectors of keypoints, and a pointer to a double array h.

* It should then compute the homography between the keypoint matches using RANSAC, store it in the array h, and return the number of inliers.

### Implementation and Result Showcase

* I have prepared a sh script to add small automation to my compile process. Therefore you can compile given files with this script. Compilation process uses valgrind to check and find any memory leak my occur and this check takes time to complete.If you want you can disable it inside of the compile.sh file.

* I have implemented methods in the problem statement and some others to help me to solve the problems more effectively.

#### Setup

* When you compile with compile.sh it runs the program once with valgrind. After that you can also run it yourself.

```bash
$ sh compile.sh
$ cd build
$ ./image-test
```

* You can play with image_test.cc to get different results. However don't forget to compile whenever you made a change :)

```C++
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
```
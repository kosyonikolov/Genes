#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <span>

#include <opencv2/opencv.hpp>

#include <genetic/opt.h>
#include <genetic/genetic.h>

struct CameraParams
{
    float focalLength;
    float ox, oy;
    std::array<float, 4> distortion;
};

std::vector<cv::Vec3f> makeGrid(const int rows, const int cols, 
                                const float step, const cv::Vec3f & offset)
{
    std::vector<cv::Vec3f> result(rows * cols);

    const float midY = static_cast<float>(rows) / 2;
    const float midX = static_cast<float>(cols) / 2;
    int i = 0;
    for (int y = 0; y < rows; y++)
    {
        const float fy = step * (y - midY);
        for (int x = 0; x < cols; x++, i++)
        {
            const float fx = step * (x - midX);
            result[i][0] = fx + offset[0];
            result[i][1] = fy + offset[1];
            result[i][2] = offset[2];
        }
    }
    
    return result;
}

void project(const std::vector<cv::Vec3f> & in, std::vector<cv::Point2f> & out, 
             const CameraParams & cam, const cv::Vec3f & rot, const cv::Vec3f & trans)
{
    float camMat0[9] = 
    {
        cam.focalLength, 0, cam.ox,
        0, cam.focalLength, cam.oy,
        0, 0, 1
    };
    cv::Mat camMat(3, 3, CV_32FC1, camMat0);

    if (out.size() != in.size())
    {
        out.resize(in.size());
    }
    cv::projectPoints(in, rot, trans, camMat, cam.distortion, out);
}

double mse(const std::span<const cv::Point2f> a, const std::span<const cv::Point2f> b)
{
    const int n = a.size();
    assert(n == b.size());
    assert(n > 0);

    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        auto dx = a[i].x - b[i].x;
        auto dy = a[i].y - b[i].y;
        sum += dx * dx + dy * dy;
    }

    return sum / n;
}

struct Solution
{
    CameraParams cam;
    cv::Vec3f trans;
    cv::Vec3f rot;
};

constexpr int xLength = 1 + 2 + 4 + 3 + 3;

void x2sol(std::span<double> x, Solution & sol)
{
    assert(x.size() == xLength);
    sol.cam.focalLength = x[0];
    sol.cam.ox = x[1];
    sol.cam.oy = x[2];
    std::copy_n(x.begin() + 3, sol.cam.distortion.size(), sol.cam.distortion.begin());
    sol.trans[0] = x[7];
    sol.trans[1] = x[8];
    sol.trans[2] = x[9];
    sol.rot[0] = x[10];
    sol.rot[1] = x[11];
    sol.rot[2] = x[12];
}

std::pair<Solution, double> solveCam(const std::vector<cv::Vec3f> & objPoints, const std::vector<cv::Point2f> & imgPoints,
                                     const int width, const int height)
{
    const auto OX_BASE = width / 2.0f;
    const auto OY_BASE = height / 2.0f;

    std::array<double, xLength> xLow  = { 600, OX_BASE - 100, OY_BASE - 100, -0.1, -0.1, -0.05, -0.05, -1, -1, -1, -0.2, -0.2, -0.2 };
    std::array<double, xLength> xHigh = { 900, OX_BASE + 100, OY_BASE + 100, +0.1, +0.1, +0.05, +0.05, +1, +1, +1, +0.2, +0.2, +0.2 };
    std::array<double, xLength> x0 = {};

    Solution tmp;
    std::vector<cv::Point2f> pts;
    auto errf = [&](const std::span<double> x)
    {
        x2sol(x, tmp);
        project(objPoints, pts, tmp.cam, tmp.rot, tmp.trans);
        auto err = mse(imgPoints, pts);
        return err;
    };

    g::Options opt;
    opt.populationSize = 20;
    opt.tournamentSize = 4;
    opt.pCrossover = 0.5;
    opt.pUniformMutation = 0.02;
    opt.big.p = 3.0 / xLength;
    opt.big.sigma = 0.1;
    opt.small.p = 0.3;
    opt.small.sigma = 0.01;

    auto gSol = g::solve(errf, x0, xLow, xHigh, opt, 500);

    std::pair<Solution, double> result;
    x2sol(gSol.x, result.first);
    result.second = std::sqrt(gSol.value);

    return result;
}

int main()
{
    const int width = 1400;
    const int height = 800;

    auto grid1 = makeGrid(11, 7, 1.5, cv::Vec3f(0.5, 0.5, 15));
    auto grid2 = makeGrid(11, 7, 0.75, cv::Vec3f(4.7, 0.3, 7));
    grid1.insert(grid1.end(), grid2.begin(), grid2.end());
    auto grid3 = makeGrid(11, 7, 0.75, cv::Vec3f(-4.2, 0.3, 7));
    grid1.insert(grid1.end(), grid3.begin(), grid3.end());

    CameraParams camTrue;
    camTrue.focalLength = 700;
    camTrue.ox = 50 + width / 2.0f;
    camTrue.oy = -15 + height / 2.0f;
    // std::fill_n(cam.distortion.begin(), cam.distortion.size(), 0);
    camTrue.distortion = { -5.0e-2f, -21.0e-3, 0.01, -0.02 };

    cv::Vec3f translationTrue(-0.25, 0.2, -0.5);
    cv::Vec3f rotationTrue(0, 0.1, 0.1);

    std::vector<cv::Point2f> imgPoints;
    project(grid1, imgPoints, camTrue, rotationTrue, translationTrue);

    auto optimSol = solveCam(grid1, imgPoints, width, height);
    std::cout << "Final err = " << optimSol.second << "\n";

    std::vector<cv::Point2f> optPoints;
    project(grid1, optPoints, optimSol.first.cam, optimSol.first.rot, optimSol.first.trans);

    cv::Mat image = cv::Mat::zeros(height, width, CV_8UC1);
    for (int i = 0; i < imgPoints.size(); i++)
    {
        cv::circle(image, imgPoints[i], 3, 255);
        cv::circle(image, optPoints[i], 3, 128);
        cv::line(image, imgPoints[i], optPoints[i], 128);
    }

    cv::imshow("Points", image);
    cv::waitKey();

    return 0;
}
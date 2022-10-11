#include <iostream>
#include "D:/DownLoad/eigen-3.4.0/Eigen/Dense"
#include <math.h>
using namespace std;
using namespace Eigen;

int main()
{
    //焦距
    double f = 0.15324; //单位m

    //已知点数据 按照x, y, X, Y, Z存储
    Array<double,1,5> point1 {-0.08615, -0.06899, 36589.41, 25273.32, 2195.17};
    Array<double,1,5> point2 {-0.05340, 0.08221, 37631.08, 31324.51, 728.69};
    Array<double,1,5> point3 {-0.01478, -0.07663, 39100.97, 24934.98, 2386.50};
    Array<double,1,5> point4 {0.01046, 0.06443, 40426.54, 30319.81, 757.31};


    //外方位元素初始值
    double Xso = (point1(0,2) + point2(0,2) + point3(0,2) + point4(0,2))/4;
    double Yso = (point1(0,3) + point2(0,3) + point3(0,3) + point4(0,3))/4;
    //航高 计算：图像中一段距离比上实地相应一段距离为 1/m 由于 1/m = f/H = l/L
    double Zso = f/(sqrt(pow(point1(0,0)-point2(0,0),2)+pow(point1(0,1)-point2(0,1),2))/
                 sqrt(pow(point1(0,2)-point2(0,2),2)+pow(point1(0,3)-point2(0,3),2)));
    //网上有人用这个
    //double Zso = (point1(0,4) + point2(0,4) + point3(0,4) + point4(0,4))/4 + f*50000;
    double fo = 0;
    double wo = 0;
    double ko = 0;

    MatrixXd X; //外方位元素改正矩阵
    MatrixXd Ro; //旋转矩阵
    while(1){//0.00000048
        /**********************************************************
         * a1,a2,a3,b1,b2,b3,c1,c2,c3为旋转矩阵的九个元素
         *********************************************************/
        //start
        double a1 = cos(fo)*cos(ko)-sin(fo)*sin(wo)*sin(ko);
        double a2 = -cos(fo)*sin(ko)-sin(fo)*sin(wo)*cos(ko);
        double a3 = -sin(fo)*cos(wo);
        double b1 = cos(wo)*sin(ko);
        double b2 = cos(wo)*cos(ko);
        double b3 = -sin(wo);
        double c1 = sin(fo)*cos(ko)+cos(fo)*sin(wo)*sin(ko);
        double c2 = -sin(fo)*sin(ko)+cos(fo)*sin(wo)*cos(ko);
        double c3 = cos(fo)*cos(wo);
        //end

        //Ro = MatrixXd {{a1,a2,a3},{b1,b2,b3},{c1,c2,c3}};

        /**********************************************************
         * 每个点对应的x，y的近似值
         *********************************************************/
        //start
        double x1js = -f*(a1*(point1(0,2)-Xso)+b1*(point1(0,3)-Yso)+c1*(point1(0,4)-Zso))/
                        (a3*(point1(0,2)-Xso)+b3*(point1(0,3)-Yso)+c3*(point1(0,4)-Zso));
        double y1js = -f*(a2*(point1(0,2)-Xso)+b2*(point1(0,3)-Yso)+c2*(point1(0,4)-Zso))/
                        (a3*(point1(0,2)-Xso)+b3*(point1(0,3)-Yso)+c3*(point1(0,4)-Zso));
        //point2
        double x2js = -f*(a1*(point2(0,2)-Xso)+b1*(point2(0,3)-Yso)+c1*(point2(0,4)-Zso))/
                        (a3*(point2(0,2)-Xso)+b3*(point2(0,3)-Yso)+c3*(point2(0,4)-Zso));
        double y2js = -f*(a2*(point2(0,2)-Xso)+b2*(point2(0,3)-Yso)+c2*(point2(0,4)-Zso))/
                        (a3*(point2(0,2)-Xso)+b3*(point2(0,3)-Yso)+c3*(point2(0,4)-Zso));
        //point3
        double x3js = -f*(a1*(point3(0,2)-Xso)+b1*(point3(0,3)-Yso)+c1*(point3(0,4)-Zso))/
                        (a3*(point3(0,2)-Xso)+b3*(point3(0,3)-Yso)+c3*(point3(0,4)-Zso));
        double y3js = -f*(a2*(point3(0,2)-Xso)+b2*(point3(0,3)-Yso)+c2*(point3(0,4)-Zso))/
                        (a3*(point3(0,2)-Xso)+b3*(point3(0,3)-Yso)+c3*(point3(0,4)-Zso));
        //point4
        double x4js = -f*(a1*(point4(0,2)-Xso)+b1*(point4(0,3)-Yso)+c1*(point4(0,4)-Zso))/
                        (a3*(point4(0,2)-Xso)+b3*(point4(0,3)-Yso)+c3*(point4(0,4)-Zso));
        double y4js = -f*(a2*(point4(0,2)-Xso)+b2*(point4(0,3)-Yso)+c2*(point4(0,4)-Zso))/
                        (a3*(point4(0,2)-Xso)+b3*(point4(0,3)-Yso)+c3*(point4(0,4)-Zso));
        //end

        /**********************************************************
         * l = x - (x)
         * L矩阵
         *********************************************************/
        //start
        double lx1 = point1(0,0) - x1js;
        double ly1 = point1(0,1) - y1js;
        double lx2 = point2(0,0) - x2js;
        double ly2 = point2(0,1) - y2js;
        double lx3 = point3(0,0) - x3js;
        double ly3 = point3(0,1) - y3js;
        double lx4 = point4(0,0) - x4js;
        double ly4 = point4(0,1) - y4js;

        Matrix<double,8,1> L;
        L << lx1,
             ly1,
             lx2,
             ly2,
             lx3,
             ly3,
             lx4,
             ly4;
        //end

        /**********************************************************
         * Z_ 为Z上一横
         * A系数矩阵
         *********************************************************/
        //start
        double Z_ = a3*(point1(0,2)-Xso)+b3*(point1(0,3)-Yso)+c3*(point1(0,4)-Zso);
        Matrix<double,2,6> A1 {{(a1*f+a3*point1(0,0))/Z_,
                                (b1*f+b3*point1(0,0))/Z_,
                                (c1*f+c3*point1(0,0))/Z_,
                                point1(0,1)*sin(wo)-(point1(0,0)/f*(point1(0,0)*cos(ko)-point1(0,1)*sin(ko))+f*cos(ko))*cos(wo),
                                -f*sin(ko)-point1(0,0)/f*(point1(0,0)*sin(ko)+point1(0,1)*cos(ko)),
                                point1(0,1)},
                               {(a2*f+a3*point1(0,1))/Z_,
                                (b2*f+b3*point1(0,1))/Z_,
                                (c2*f+c3*point1(0,1))/Z_,
                                -point1(0,0)*sin(wo)-(point1(0,1)/f*(point1(0,0)*cos(ko)-point1(0,1)*sin(ko))-f*sin(ko))*cos(wo),
                                -f*cos(ko)-point1(0,1)/f*(point1(0,0)*sin(ko)+point1(0,1)*cos(ko)),
                                -point1(0,0)}};

        Matrix<double,2,6> A2 {{(a1*f+a3*point2(0,0))/Z_,
                                (b1*f+b3*point2(0,0))/Z_,
                                (c1*f+c3*point2(0,0))/Z_,
                                point2(0,1)*sin(wo)-(point2(0,0)/f*(point2(0,0)*cos(ko)-point2(0,1)*sin(ko))+f*cos(ko))*cos(wo),
                                -f*sin(ko)-point2(0,0)/f*(point2(0,0)*sin(ko)+point2(0,1)*cos(ko)),
                                point2(0,1)},
                               {(a2*f+a3*point2(0,1))/Z_,
                                (b2*f+b3*point2(0,1))/Z_,
                                (c2*f+c3*point2(0,1))/Z_,
                                -point2(0,0)*sin(wo)-(point2(0,1)/f*(point2(0,0)*cos(ko)-point2(0,1)*sin(ko))-f*sin(ko))*cos(wo),
                                -f*cos(ko)-point2(0,1)/f*(point2(0,0)*sin(ko)+point2(0,1)*cos(ko)),
                                -point2(0,0)}};
        Matrix<double,2,6> A3 {{(a1*f+a3*point3(0,0))/Z_,
                                (b1*f+b3*point3(0,0))/Z_,
                                (c1*f+c3*point3(0,0))/Z_,
                                point3(0,1)*sin(wo)-(point3(0,0)/f*(point3(0,0)*cos(ko)-point3(0,1)*sin(ko))+f*cos(ko))*cos(wo),
                                -f*sin(ko)-point3(0,0)/f*(point3(0,0)*sin(ko)+point3(0,1)*cos(ko)),
                                point3(0,1)},
                               {(a2*f+a3*point3(0,1))/Z_,
                                (b2*f+b3*point3(0,1))/Z_,
                                (c2*f+c3*point3(0,1))/Z_,
                                -point3(0,0)*sin(wo)-(point3(0,1)/f*(point3(0,0)*cos(ko)-point3(0,1)*sin(ko))-f*sin(ko))*cos(wo),
                                -f*cos(ko)-point3(0,1)/f*(point3(0,0)*sin(ko)+point3(0,1)*cos(ko)),
                                -point3(0,0)}};
        Matrix<double,2,6> A4 {{(a1*f+a3*point4(0,0))/Z_,
                                (b1*f+b3*point4(0,0))/Z_,
                                (c1*f+c3*point4(0,0))/Z_,
                                point4(0,1)*sin(wo)-(point4(0,0)/f*(point4(0,0)*cos(ko)-point4(0,1)*sin(ko))+f*cos(ko))*cos(wo),
                                -f*sin(ko)-point4(0,0)/f*(point4(0,0)*sin(ko)+point4(0,1)*cos(ko)),
                                point4(0,1)},
                               {(a2*f+a3*point4(0,1))/Z_,
                                (b2*f+b3*point4(0,1))/Z_,
                                (c2*f+c3*point4(0,1))/Z_,
                                -point4(0,0)*sin(wo)-(point4(0,1)/f*(point4(0,0)*cos(ko)-point4(0,1)*sin(ko))-f*sin(ko))*cos(wo),
                                -f*cos(ko)-point4(0,1)/f*(point4(0,0)*sin(ko)+point4(0,1)*cos(ko)),
                                -point4(0,0)}};
        Matrix<double,8,6> A;
        A << A1,
             A2,
             A3,
             A4;
        //end


        /**********************************************************
         * X为外方位元素改正数
         *********************************************************/
        X = (A.transpose()*A).inverse()*A.transpose()*L;
        //V = A*X-L;

        Xso += X(0,0);
        Yso += X(1,0);
        Zso += X(2,0);
        fo += X(3,0);
        wo += X(4,0);
        ko += X(5,0);

        if(abs(X(5,0)) < 0.00000048 && abs(X(4,0)) < 0.00000048 && abs(X(3,0)) < 0.00000048){
            break;
        }

    }

    cout << "The outer bearing element" << endl;
    cout << "Xso:" << Xso << " " << "Yso:" << Yso << " " << "Zso:" << Zso << endl << "fo:" << fo << " " << "wo:" << wo << " " << "ko:" << ko << endl;
    return 0;
}

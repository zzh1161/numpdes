#include<fstream>
#include"../include/json.hpp"
#include"../include/GenericFactory.hpp"
#include"../include/IncAllMethod.h"

using TimeIntegratorFactory = GenericFactory<TimeIntegrator<6>, std::string, 
									std::unique_ptr<TimeIntegrator<6>>(*)(int)>;

using vec6Func = std::function<Vec<real,6>(Vec<real,6>,real)>;
using json  = nlohmann::json;

void registerFactory()
{
	TimeIntegratorFactory &Tf = TimeIntegratorFactory::createFactory();
    Tf.registerProduct("ForwardEuler",    createForwardEuler);
	Tf.registerProduct("Adams-Bashforth", createAdamsBashf);
	Tf.registerProduct("Adams-Moulton",   createAdamsMoulton);
	Tf.registerProduct("ClassicalRK",     createClassiRK);
	Tf.registerProduct("BDF",             createBDF);
	Tf.registerProduct("ESDIRK",          createESDIRK);
	Tf.registerProduct("Gauss-Legendre",  createGaussLege);
	Tf.registerProduct("F45RK",           createF45RK);
	Tf.registerProduct("DP54RK",          createDP54RK);
}

int main()
{
    json j;
    std::ifstream jfile("./src/inputfile/condition.json");
    jfile >> j;
    jfile.close();

    std::ofstream outfile("./output/plot.txt", std::ios::out);

    Vec<real,6> U0_1{0.994, 0, 0, 0, -2.0015851063790825224, 0};
    real T1 = 17.06521656015796;
    Vec<real,6> U0_2{0.87978, 0, 0, 0, -0.3797, 0};
    real T2 = 19.14045706162071;
    vec6Func f = [](Vec<real,6> u, real t){
        real mu = 1.0/81.45;
        real f1 = u[3];
        real f2 = u[4];
        real f3 = u[5];
        real f4 = 2*u[4]+u[0] - (mu*(u[0]+mu-1))/pow((pow(u[1],2)+pow(u[2],2)+pow(u[0]+mu-1,2)),1.5)
                  - ((1-mu)*(u[0]+mu))/pow((pow(u[1],2)+pow(u[2],2)+pow(u[0]+mu,2)),1.5);
        real f5 = -2*u[3]+u[1] - (mu*u[1])/pow((pow(u[1],2)+pow(u[2],2)+pow(u[0]+mu-1,2)),1.5)
                  - ((1-mu)*u[1])/pow((pow(u[1],2)+pow(u[2],2)+pow(u[0]+mu,2)),1.5);
        real f6 = -(mu*u[2])/pow((pow(u[1],2)+pow(u[2],2)+pow(u[0]+mu-1,2)),1.5)
                  -((1-mu)*u[2])/pow((pow(u[1],2)+pow(u[2],2)+pow(u[0]+mu,2)),1.5);
        return Vec<real,6>{f1,f2,f3,f4,f5,f6};
    };

	registerFactory();
    TimeIntegratorFactory &tf = TimeIntegratorFactory::createFactory();

    clock_t begin, end;

    outfile << "!!! Forward Euler !!!" << std::endl;
    auto FacForEuler = tf.createProduct("ForwardEuler", j["FE"]["steps"]);
    outfile << "(10.198) with T = T1" << std::endl;
    begin = clock();
    FacForEuler->solve(U0_1, f, T1);
    end = clock();
    std::cout << "Euler method for (10.198) takes " << 
                (end-begin)*1.0/CLOCKS_PER_SEC << " seconds." << std::endl;
    outfile << *FacForEuler << std::endl;
    outfile << "(10.199) with T = T2" << std::endl;
    begin = clock();
    FacForEuler->solve(U0_2, f, T2);
    end = clock();
    std::cout << "Euler method for (10.199) takes " << 
                (end-begin)*1.0/CLOCKS_PER_SEC << " seconds." << std::endl;
    outfile << *FacForEuler << std::endl;
    outfile << std::endl;

    outfile << "!!! Classical Runge-Kutta !!!" << std::endl;
    auto FacClassRK = tf.createProduct("ClassicalRK", j["cRK"]["steps"]);
    outfile << "(10.198) with T = T1" << std::endl;
    begin = clock();
    FacClassRK->solve(U0_1, f, T1);
    end = clock();
    std::cout << "Classical RK method for (10.198) takes " << 
                (end-begin)*1.0/CLOCKS_PER_SEC << " seconds." << std::endl;
    outfile << *FacClassRK << std::endl;
    outfile << "(10.199) with T = T2" << std::endl;
    begin = clock();
    FacClassRK->solve(U0_2, f, T2);
    end = clock();
    std::cout << "Classical RK method for (10.199) takes " << 
                (end-begin)*1.0/CLOCKS_PER_SEC << " seconds." << std::endl;
    outfile << *FacClassRK << std::endl;
    outfile << std::endl;

    outfile << "!!! Dormand-Prince 5(4) embedded RK !!!" << std::endl;
    auto FacDP54RK = tf.createProduct("DP54RK", j["DP54"]["steps"]);
    outfile << "(10.198) with T = T1" << std::endl;
    begin = clock();
    FacDP54RK->solve(U0_1, f, T1);
    end = clock();
    std::cout << "Dormand-Prince 5(4) RK method for (10.198) takes " << 
                (end-begin)*1.0/CLOCKS_PER_SEC << " seconds." << std::endl;
    outfile << *FacDP54RK << std::endl;
    outfile << "(10.199) with T = T2" << std::endl;
    begin = clock();
    FacDP54RK->solve(U0_2, f, T2);
    end = clock();
    std::cout << "Dormand-Prince 5(4) RK method for (10.199) takes " << 
                (end-begin)*1.0/CLOCKS_PER_SEC << " seconds." << std::endl;
    outfile << *FacDP54RK << std::endl;
    outfile << std::endl;

	return 0;
}
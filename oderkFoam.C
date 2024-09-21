#include "fvCFD.H"
#include "pisoControl.H"
#include "fvOptions.H"

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createPhi.H"

    Info<< "\nStarting time loop\n" << endl;

    dimensionedScalar b("b", dimensionSet(0, 1, -1, 0, 0, 0, 0), 1.0);
    double a = 1.0;
    double omega = 10.0;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        volVectorField Uold = U;
        
        scalar t = runTime.value(); // current time 
        scalar dt = runTime.deltaT().value(); // time step size

        volVectorField k1 = dt*(-a*U + b*vector(Foam::cos(omega*t), 0, 0));
        volVectorField k2 = dt*(-a*(U + 0.5*k1) + b*vector(Foam::cos(omega*(t + 0.5*dt)), 0, 0));
        volVectorField k3 = dt*(-a*(U + 0.5*k2) + b*vector(Foam::cos(omega*(t + 0.5*dt)), 0, 0));
        volVectorField k4 = dt*(-a*(U + k3) + b*vector(Foam::cos(omega*(t + dt)), 0, 0));

        U = U + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);

        U.correctBoundaryConditions();

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}
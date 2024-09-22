#include "fvCFD.H"

// ODE function
doubleScalar f(const doubleScalar& t, const doubleScalar& y, const doubleScalar& a, const doubleScalar& b, const doubleScalar& omega)
{
    Info << "y in the ode: " << y << nl;
    return -a * y + b * Foam::cos(omega * t);
}

// Runge-Kutta 4th order step
doubleScalar runge_kutta_step(const doubleScalar& t, const doubleScalar& y, const doubleScalar& h,
                              const doubleScalar& a, const doubleScalar& b, const doubleScalar& omega)
{
    doubleScalar k1 = h * f(t, y, a, b, omega);
    Info << " t: " << t << nl;
    Info << "ODE K1: " << f(t, y, a, b, omega) << nl;
    Info << "k1: " << k1 << nl;

    doubleScalar k2 = h * f(t + 0.5*h, y + 0.5*k1, a, b, omega);
    Info << "k2: " << k2 << nl;
    Info << "ODE K2 " << f(t + 0.5*h, y + 0.5*k1, a, b, omega) << nl;

    doubleScalar k3 = h * f(t + 0.5*h, y + 0.5*k2, a, b, omega);
    Info << "k3: " << k3 << nl;
    Info << "ODE K3 " << f(t + 0.5*h, y + 0.5*k2, a, b, omega) << nl;

    doubleScalar k4 = h * f(t + h, y + k3, a, b, omega);
    Info << "k4: " << k4 << nl;
    Info << "ODE K4 " << f(t + h, y + k3, a, b, omega) << nl;

    Info << "y: " << y << nl;

    return y + (k1 + 2*k2 + 2*k3 + k4) / 6;
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Read transport properties
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Read a, b, and omega from transportProperties
    dimensionedScalar a("a", dimless, transportProperties);
    dimensionedScalar b("b", dimVelocity, transportProperties);
    dimensionedScalar omega("omega", dimless/dimTime, transportProperties);

    // Create the U field
    volScalarField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Set initial condition
    scalar y0 = a.value() * b.value() / (a.value()*a.value() + omega.value()*omega.value());
    U.primitiveFieldRef() = y0;

    // Main time loop
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Get current time and time step
        doubleScalar t = runTime.value() - runTime.deltaTValue(); // Use previous time step
        doubleScalar dt = runTime.deltaTValue();

        forAll(U, cellI)
        {
            doubleScalar y = U[cellI];
            U[cellI] = runge_kutta_step(t, y, dt, a.value(), b.value(), omega.value());
            Info << "y_rk[" << runTime.timeIndex() << "]: " << U[cellI] << nl;
            Info << "**********" << nl;
        }

        U.correctBoundaryConditions();

        runTime.write();

        Info<< "Time: " << runTime.timeName() << ", U: " << U[0] << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
package Integrals;

import Molecule.Atom;
import Molecule.MoleculeIntegral;

import java.util.List;

public class OneIntegrals {
    private final double[][] overlapInt;
    private final double[][] kineticInt;

    // Constructor — builds overlap integral matrix automatically
    public OneIntegrals(MoleculeIntegral molInt) {
        this.overlapInt = computeOverlap(molInt);
        this.kineticInt = computeKinetic(molInt);
    }

    // Getter integrals matrix
    public double[][] getOverlapInt() {
        return overlapInt;
    }
    public double[][] getKineticInt() {
        return kineticInt;
    }

    public static double[][] computeOverlap(MoleculeIntegral molInt) {

        List<Atom.AtomIntegrals> atomInt = molInt.getAtoms();
        int nOrb = atomInt.size();

        double[][] intOvelap = new double[nOrb][nOrb];

        for (int i = 0; i < nOrb; i++) {
            for (int j = i; j < nOrb; j++) {
                intOvelap[i][j] = overlap(atomInt.get(i),atomInt.get(j));
                intOvelap[j][i] = intOvelap[i][j];
            }
        }
        return intOvelap;
    }

    public static double[][] computeKinetic(MoleculeIntegral molInt) {

        List<Atom.AtomIntegrals> atomInt = molInt.getAtoms();
        int nOrb = atomInt.size();

        double[][] intKinetic = new double[nOrb][nOrb];

        for (int i = 0; i < nOrb; i++) {
            for (int j = i; j < nOrb; j++) {
                intKinetic[i][j] = kinetic(atomInt.get(i),atomInt.get(j));
                intKinetic[j][i] = intKinetic[i][j];
            }
        }
        return intKinetic;
    }

    public static double overlap(Atom.AtomIntegrals atomInt1,
                                 Atom.AtomIntegrals atomInt2) {

        double resultado = 0.0;

        for (int i = 0; i < atomInt1.exponent().size(); i++) {
            for (int j = 0; j < atomInt2.exponent().size(); j++) {
                resultado += normOrbital(atomInt1.exponent().get(i), atomInt1.cartAngular())
                        * normOrbital(atomInt2.exponent().get(j), atomInt2.cartAngular())
                        * atomInt1.coefficient().get(i) * atomInt2.coefficient().get(j)
                        * ovInt(atomInt1.exponent().get(i), atomInt2.exponent().get(j), atomInt1.coord(),
                        atomInt2.coord(), atomInt1.cartAngular(), atomInt2.cartAngular());
            }
        }

        return resultado;
    }

    public static double kinetic(Atom.AtomIntegrals atomInt1,
                                 Atom.AtomIntegrals atomInt2) {

        double resultado = 0.0;

        for (int i = 0; i < atomInt1.exponent().size(); i++) {
            for (int j = 0; j < atomInt2.exponent().size(); j++) {
                resultado += normOrbital(atomInt1.exponent().get(i), atomInt1.cartAngular())
                        * normOrbital(atomInt2.exponent().get(j), atomInt2.cartAngular())
                        * atomInt1.coefficient().get(i) * atomInt2.coefficient().get(j)
                        * kinInt(atomInt1.exponent().get(i), atomInt2.exponent().get(j), atomInt1.coord(),
                        atomInt2.coord(), atomInt1.cartAngular(), atomInt2.cartAngular());
            }
        }

        return resultado;
    }

    public static double ovInt(double orbC1, double orbC2,
                               Atom.Coordinates coord1, Atom.Coordinates coord2,
                               Atom.AngCoordinates angCoord1, Atom.AngCoordinates angCoord2) {

        double producto = 1.0;

        // Product over x, y, z components
        producto *= sRecurr(orbC1, orbC2, coord1.x(), coord2.x(), angCoord1.x(), angCoord2.x());
        producto *= sRecurr(orbC1, orbC2, coord1.y(), coord2.y(), angCoord1.y(), angCoord2.y());
        producto *= sRecurr(orbC1, orbC2, coord1.z(), coord2.z(), angCoord1.z(), angCoord2.z());

        // Compute distance squared between coord1 and coord2
        double dist2 = Math.pow(coord1.x() - coord2.x(), 2) +
                Math.pow(coord1.y() - coord2.y(), 2) +
                Math.pow(coord1.z() - coord2.z(), 2);

        double eAB = Math.exp(-(orbC1 * orbC2 / (orbC1 + orbC2)) * dist2);

        return eAB * Math.pow(Math.PI / (orbC1 + orbC2), 1.5) * producto;
    }

    public static double kinInt(double orbC1, double orbC2,
                               Atom.Coordinates coord1, Atom.Coordinates coord2,
                               Atom.AngCoordinates angCoord1, Atom.AngCoordinates angCoord2) {


        double producto1 = kiRecurr(orbC1, orbC2, coord1.x(), coord2.x(), angCoord1.x(), angCoord2.x())
                * sRecurr(orbC1, orbC2, coord1.y(), coord2.y(), angCoord1.y(), angCoord2.y())
                * sRecurr(orbC1, orbC2, coord1.z(), coord2.z(), angCoord1.z(), angCoord2.z());

        double producto2 = sRecurr(orbC1, orbC2, coord1.x(), coord2.x(), angCoord1.x(), angCoord2.x())
                * kiRecurr(orbC1, orbC2, coord1.y(), coord2.y(), angCoord1.y(), angCoord2.y())
                * sRecurr(orbC1, orbC2, coord1.z(), coord2.z(), angCoord1.z(), angCoord2.z());

        double producto3 = sRecurr(orbC1, orbC2, coord1.x(), coord2.x(), angCoord1.x(), angCoord2.x())
                * sRecurr(orbC1, orbC2, coord1.y(), coord2.y(), angCoord1.y(), angCoord2.y())
                * kiRecurr(orbC1, orbC2, coord1.z(), coord2.z(), angCoord1.z(), angCoord2.z());

        double suma = producto1 + producto2 + producto3;

        // Distance squared between centers
        double dist2 = Math.pow(coord1.x() - coord2.x(), 2) +
                Math.pow(coord1.y() - coord2.y(), 2) +
                Math.pow(coord1.z() - coord2.z(), 2);

        double EAB = Math.exp(-(orbC1 * orbC2 / (orbC1 + orbC2)) * dist2);

        return EAB * Math.pow(Math.PI / (orbC1 + orbC2), 1.5) * suma;

    }

    public static double sRecurr(double alpha1, double beta1,
                                 double coord1, double coord2,
                                 int angCoord1coord, int angCoord2coord) {

        // Base cases
        if (angCoord1coord == 0 && angCoord2coord == 0) {
            return 1.0;
        } else if (angCoord1coord == 1 && angCoord2coord == 0) {
            return -(coord1 - ((alpha1 * coord1 + beta1 * coord2) / (alpha1 + beta1)));
        } else if (angCoord1coord > 1 && angCoord2coord == 0) {
            double term1 = -(coord1 - ((alpha1 * coord1 + beta1 * coord2) / (alpha1 + beta1)))
                    * sRecurr(alpha1, beta1, coord1, coord2, angCoord1coord - 1, angCoord2coord);
            double term2 = ((angCoord1coord - 1) / (2.0 * (alpha1 + beta1)))
                    * sRecurr(alpha1, beta1, coord1, coord2, angCoord1coord - 2, angCoord2coord);
            return term1 + term2;
        } else {
            return sRecurr(alpha1, beta1, coord1, coord2, angCoord1coord + 1, angCoord2coord - 1)
                    + (coord1 - coord2)
                    * sRecurr(alpha1, beta1, coord1, coord2, angCoord1coord, angCoord2coord - 1);
        }
    }

    private static double kiRecurr(double alpha, double beta,
                                        double coord1, double coord2,
                                        int CA1, int CA2) {

        double resultado = 0.0;

        // Case 1: both angular momenta are zero
        if (CA1 == 0 && CA2 == 0) {
            resultado = 2.0 * alpha * beta * sRecurr(alpha, beta, coord1, coord2, 1, 1);

            // Case 2: CA1 > 0, CA2 == 0
        } else if (CA1 > 0 && CA2 == 0) {
            resultado = -CA1 * beta * sRecurr(alpha, beta, coord1, coord2, CA1 - 1, 1)
                    + 2.0 * alpha * beta * sRecurr(alpha, beta, coord1, coord2, CA1 + 1, 1);

            // Case 3: CA1 == 0, CA2 > 0
        } else if (CA1 == 0 && CA2 > 0) {
            resultado = -CA2 * alpha * sRecurr(alpha, beta, coord1, coord2, 1, CA2 - 1)
                    + 2.0 * alpha * beta * sRecurr(alpha, beta, coord1, coord2, 1, CA2 + 1);

            // Case 4: both CA1 and CA2 > 0
        } else {
            resultado = (
                    (CA1 * CA2 * sRecurr(alpha, beta, coord1, coord2, CA1 - 1, CA2 - 1))
                            - (2.0 * CA1 * beta * sRecurr(alpha, beta, coord1, coord2, CA1 - 1, CA2 + 1))
                            - (2.0 * CA2 * alpha * sRecurr(alpha, beta, coord1, coord2, CA1 + 1, CA2 - 1))
                            + (4.0 * alpha * beta * sRecurr(alpha, beta, coord1, coord2, CA1 + 1, CA2 + 1))
            ) / 2.0;
        }

        return resultado;
    }

    public static double normOrbital(double orbC, Atom.AngCoordinates vecA) {

        double pi = Math.PI;
        double b = Math.pow((2.0 * orbC / pi), 0.75);
        b *= Math.pow((4 * orbC), (vecA.x() + vecA.y() + vecA.z()) / 2.0);

        b /= Math.sqrt(
                doubleFactorial(2 * vecA.x() - 1)
                        * doubleFactorial(2 * vecA.y() - 1)
                        * doubleFactorial(2 * vecA.z() - 1)
        );

        return b;
    }

    public static double doubleFactorial(int n) {
        if (n <= 0) return 1.0;
        double result = 1.0;
        for (int i = n; i > 0; i -= 2) {
            result *= i;
        }
        return result;
    }
}

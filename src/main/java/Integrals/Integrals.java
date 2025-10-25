package Integrals;

import Molecule.Atom;
import Molecule.MoleculeIntegral;

import java.util.List;

public class Integrals {

    public static double[][] intSolapamiento(MoleculeIntegral molInt) {


        List<Atom.AtomIntegrals> atomInt = molInt.getAtoms();
        int nOrb = atomInt.size();

        double[][] intSolap = new double[nOrb][nOrb];
        for (Atom.AtomIntegrals atom : molInt.getAtoms()) {
            System.out.printf("%-3s %15.6f %15.6f %15.6f%n",
                    atom.symbol(), atom.coord().x(), atom.coord().y(), atom.coord().z());
            System.out.printf("     Orbital: %-4s%n", atom.orbital());
            System.out.printf("     Exponents:   %s%n", atom.exponent());
            System.out.printf("     Coefficients:%s%n", atom.coefficient());
            System.out.printf("     Cart. Ang.:  %s%n%n", atom.cartAngular());
        }

        for (int i = 0; i < nOrb; i++) {
            for (int j = i; j < nOrb; j++) {
                intSolap[i][j] = solapamiento(atomInt.get(i),atomInt.get(j));
                intSolap[j][i] = intSolap[i][j];
            }
        }
        return intSolap;
    }

    public static double solapamiento(Atom.AtomIntegrals atomInt1,
                                      Atom.AtomIntegrals atomInt2) {

        double resultado = 0.0;

        // Double loop: i=1,a_fin ; j=1,b_fin  (Fortran indices start at 1)
        for (int i = 0; i < atomInt1.exponent().size(); i++) {
            for (int j = 0; j < atomInt2.exponent().size(); j++) {
                resultado += normOrbital(atomInt1.coefficient().get(i), atomInt1.cartAngular())
                        * normOrbital(atomInt2.coefficient().get(j), atomInt2.cartAngular())
                        * atomInt1.exponent().get(i) * atomInt2.exponent().get(j)
                        * Ov_int(atomInt1.coefficient().get(i), atomInt2.coefficient().get(j), atomInt1.coord(),
                        atomInt2.coord(), atomInt1.cartAngular(), atomInt2.cartAngular());
            }
        }

        return resultado;
    }

    // === Function: Ov_int ===
    public static double Ov_int(double orb_c_1, double orb_c_2,
                                Atom.coordinates coord1, Atom.coordinates coord2,
                                Atom.angCoordinates CA1, Atom.angCoordinates CA2) {

        double producto = 1.0;

        // Product over x, y, z components
        producto *= s_re(orb_c_1, orb_c_2, coord1.x(), coord2.x(), CA1.x(), CA2.x());
        producto *= s_re(orb_c_1, orb_c_2, coord1.y(), coord2.y(), CA1.y(), CA2.y());
        producto *= s_re(orb_c_1, orb_c_2, coord1.z(), coord2.z(), CA1.z(), CA2.z());

        // Compute distance squared between coord1 and coord2
        double dist2 = Math.pow(coord1.x() - coord2.x(), 2) +
                Math.pow(coord1.y() - coord2.y(), 2) +
                Math.pow(coord1.z() - coord2.z(), 2);

        double EAB = Math.exp(-(orb_c_1 * orb_c_2 / (orb_c_1 + orb_c_2)) * dist2);
        double overlap = EAB * Math.pow(Math.PI / (orb_c_1 + orb_c_2), 1.5) * producto;

        return overlap;
    }

    public static double s_re(double alpha1, double beta1,
                              double coord1_1, double coord1_2,
                              int CA1_1, int CA1_2) {

        // Base cases
        if (CA1_1 == 0 && CA1_2 == 0) {
            return 1.0;
        } else if (CA1_1 == 1 && CA1_2 == 0) {
            return -(coord1_1 - ((alpha1 * coord1_1 + beta1 * coord1_2) / (alpha1 + beta1)));
        } else if (CA1_1 > 1 && CA1_2 == 0) {
            double term1 = -(coord1_1 - ((alpha1 * coord1_1 + beta1 * coord1_2) / (alpha1 + beta1)))
                    * s_re(alpha1, beta1, coord1_1, coord1_2, CA1_1 - 1, CA1_2);
            double term2 = ((CA1_1 - 1) / (2.0 * (alpha1 + beta1)))
                    * s_re(alpha1, beta1, coord1_1, coord1_2, CA1_1 - 2, CA1_2);
            return term1 + term2;
        } else {
            return s_re(alpha1, beta1, coord1_1, coord1_2, CA1_1 + 1, CA1_2 - 1)
                    + (coord1_1 - coord1_2)
                    * s_re(alpha1, beta1, coord1_1, coord1_2, CA1_1, CA1_2 - 1);
        }
    }

    public static double normOrbital(double orbC, Atom.angCoordinates vecA) {

        double pi = Math.PI;
        double b = Math.pow((2.0 * orbC / pi), 0.75);
        b *= Math.pow((4 * orbC), (vecA.x() + vecA.y() + vecA.z()) / 2.0);

        b /= Math.sqrt(
                doubleFactorial(2 * vecA.x() - 1)
                        * doubleFactorial(2 * vecA.y() - 1)
                        * doubleFactorial(2 * vecA.z() - 1)
        );

//        boolean condition =
//                (vecA.get(0) > 0 && vecA.get(1) > 0 && vecA.get(2) == 0)
//                        || (vecA.get(0) > 0 && vecA.get(1) == 0 && vecA.get(2) > 0)
//                        || (vecA.get(0) == 0 && vecA.get(1) > 0 && vecA.get(2) > 0)
//                        || (vecA.get(0) > 0 && vecA.get(1) > 0 && vecA.get(2) > 0);
//
//        if (condition) {
//            b /= Math.sqrt(
//                    doubleFactorial(2 * (vecA.get(0) + vecA.get(1) + vecA.get(2)) - 1)
//                            / (doubleFactorial(2 * vecA.get(0) - 1)
//                            * doubleFactorial(2 * vecA.get(1) - 1)
//                            * doubleFactorial(2 * vecA.get(2) - 1))
//            );
//        }

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

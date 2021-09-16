use crate::arith::misc::*;
use crate::arith::modular_arithmetic::*;
use rug::Integer;

#[derive(Clone, PartialEq, Debug, Default)]
pub struct MontgomeryPoint {
    pub x: Integer,
    pub z: Integer,
    a24: Integer,
    modulo: Integer,
}

impl MontgomeryPoint {
    /// Montgomery Point
    #[allow(dead_code)]
    pub fn new(x: Integer, z: Integer, a24: Integer, modulo: Integer) -> MontgomeryPoint {
        MontgomeryPoint { x, z, a24, modulo }
    }

    pub fn new2(x: Integer, z: Integer, a: Integer, modulo: Integer) -> MontgomeryPoint {
        let inv = invert_mod(&Integer::from(4), &modulo).unwrap();
        MontgomeryPoint {
            x,
            z,
            a24: multiply_mod(&Integer::from(&a + 2), &inv, &modulo),
            modulo,
        }
    }

    /// Two points are equal if their ratio x.z^{-1} are congruent mod n
    #[allow(dead_code)]
    pub fn equals(&self, other: &MontgomeryPoint) -> bool {
        if !self.modulo.eq(&other.modulo) || !self.a24.eq(&other.a24) {
            false;
        }
        // Compute the inverse of z mod n...
        let self_z_inverse = invert_mod(&self.z, &self.modulo);
        let other_z_inverse = invert_mod(&other.z, &other.modulo);
        // ... provided it exists:
        match (self_z_inverse, other_z_inverse) {
            (Some(self_z_inv), Some(other_z_inv)) => {
                let self_ratio = multiply_mod(&self.x, &self_z_inv, &self.modulo);
                let other_ratio = multiply_mod(&other.x, &other_z_inv, &self.modulo);
                self_ratio == other_ratio // compare ratios x.z^{-1} mod n
            }
            _ => false, // z isn't invertible mod n
        }
    }

    pub fn addh(&self, other: &MontgomeryPoint, diff: &MontgomeryPoint) -> MontgomeryPoint {
        // diff = self - other
        // TODO: Check a24 and modulo
        let self_x_min_z = subtract_mod(&self.x, &self.z, &self.modulo);
        let self_x_plus_z = add_mod(&self.x, &self.z, &self.modulo);

        let other_x_min_z = subtract_mod(&other.x, &other.z, &self.modulo);
        let other_x_plus_z = add_mod(&other.x, &other.z, &self.modulo);

        let prod1 = multiply_mod(&self_x_min_z, &other_x_plus_z, &self.modulo);
        let prod2 = multiply_mod(&self_x_plus_z, &other_x_min_z, &self.modulo);

        let addition = add_mod(&prod1, &prod2, &self.modulo);
        let subtraction = subtract_mod(&prod1, &prod2, &self.modulo);

        let sqr1 = multiply_mod(&addition, &addition, &self.modulo);
        let sqr2 = multiply_mod(&subtraction, &subtraction, &self.modulo);

        MontgomeryPoint {
            x: multiply_mod(&diff.z, &sqr1, &self.modulo),
            z: multiply_mod(&diff.x, &sqr2, &self.modulo),
            a24: self.a24.clone(),
            modulo: self.modulo.clone(),
        }
    }

    /// Doubles a point in Montgomery form, requires five multiplications
    pub fn double(&self) -> MontgomeryPoint {
        let self_x_plus_z = add_mod(&self.x, &self.z, &self.modulo);
        let self_x_min_z = subtract_mod(&self.x, &self.z, &self.modulo);

        let u = self_x_plus_z.square();
        let v = self_x_min_z.square();

        let diff = Integer::from(&u - &v);
        let x = multiply_mod(&u, &v, &self.modulo);
        let z = take_mod(
            &Integer::from(&diff * Integer::from(&v + &self.a24 * &diff)),
            &self.modulo,
        );

        MontgomeryPoint {
            x,
            z,
            a24: Integer::from(&self.a24),
            modulo: Integer::from(&self.modulo),
        }
    }

    /// Scalar multiplication in Montgomery form
    pub fn montgomery_ladder(&self, k: &Integer) -> MontgomeryPoint {
        let mut q = self.clone();
        let mut p = self.double();
        let bv = bits(k);
        for b in 1..bv.len() {
            if bv[b] == '1' {
                q = p.addh(&q, &self);
                p = p.double();
            } else {
                p = q.addh(&p, &self);
                q = q.double();
            }
        }
        q
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn montgomery_addh_tests() {
        let p1 = MontgomeryPoint::new(
            Integer::from(11),
            Integer::from(16),
            Integer::from(7),
            Integer::from(29),
        );
        let p2 = MontgomeryPoint::new(
            Integer::from(13),
            Integer::from(10),
            Integer::from(7),
            Integer::from(29),
        );
        let p3 = p2.addh(&p1, &p1);
        assert_eq!((p3.x, p3.z), (Integer::from(23), Integer::from(17)));
    }

    #[test]
    fn montgomery_double_tests() {
        let p = MontgomeryPoint::new(
            Integer::from(11),
            Integer::from(16),
            Integer::from(7),
            Integer::from(29),
        );
        let q = p.double();
        assert_eq!((q.x, q.z), (Integer::from(13), Integer::from(10)));

        let p1 = MontgomeryPoint::new2(
            Integer::from(10),
            Integer::from(17),
            Integer::from(10),
            Integer::from(101),
        );
        let p2 = p1.double();
        assert_eq!(
            p2,
            MontgomeryPoint::new2(
                Integer::from(68),
                Integer::from(56),
                Integer::from(10),
                Integer::from(101)
            )
        );
    }

    #[test]
    fn montgomery_ladder_tests() {
        let p = MontgomeryPoint::new(
            Integer::from(11),
            Integer::from(16),
            Integer::from(7),
            Integer::from(29),
        );
        let q = p.montgomery_ladder(&Integer::from(3));
        assert_eq!((q.x, q.z), (Integer::from(23), Integer::from(17)));
    }

    #[test]
    fn montgomery_double_tests2() {
        let x = Integer::from(10);
        let z = Integer::from(17);
        let a = Integer::from(10);
        let modulo = Integer::from(101);
        let a24 = multiply_mod(
            &add_mod(&a, &Integer::from(2), &modulo),
            &invert_mod(&Integer::from(4), &modulo).unwrap(),
            &modulo,
        );
        let a24_1 = Integer::from(&a24);
        let a24_2 = Integer::from(&a24);
        let p1 = MontgomeryPoint {
            x,
            z,
            a24: a24_1,
            modulo,
        };

        let mod_2 = Integer::from(101);
        let x1 = Integer::from(68);
        let z1 = Integer::from(56);
        let p2 = p1.double();

        assert_eq!(
            p2,
            MontgomeryPoint {
                x: x1,
                z: z1,
                a24: a24_2,
                modulo: mod_2
            }
        );
    }
}

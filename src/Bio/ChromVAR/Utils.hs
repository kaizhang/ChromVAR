{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables   #-}
module Bio.ChromVAR.Utils
    ( weight
    , Whitening(..)
    , whiten
    ) where

import qualified Eigen.Matrix as E
import Eigen.Matrix (Matrix)
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Storable as S
import Statistics.Sample (mean)
import Eigen.Solver.LA
import Data.Singletons

weight :: (Double, Double) -> (Double, Double) -> Double
weight (x1, y1) (x2, y2) = gaussian $ sqrt $ (x1 - x2)**2 + (y1 - y2)**2
  where
    gaussian d = let sigma = 0.1 in exp $ negate $ d**2 / (2 * sigma**2)
    --gaussian d = d
{-# INLINE weight #-}

data Whitening = ZCA
               | Cholesky

-- | Rows are samples, columns are features.
whiten :: (SingI n, SingI m) => Whitening -> Matrix n m Double -> Matrix n m Double
whiten method mat = case method of
    Cholesky -> E.transpose $ E.inverse (cholesky cov) `E.mul` E.transpose mat
  where
    cov = covariance mat

covariance :: forall n m. (SingI n, SingI m) => Matrix n m Double -> Matrix m m Double
covariance mat = E.map (/n) $ E.transpose cs `E.mul` cs
  where
    n = fromIntegral $ E.rows mat - 1
    cs = E.fromColumns $ map f $ E.toColumns mat :: Matrix n m Double
    f x = let m = mean x in S.map (subtract m) x
{-# INLINE covariance #-}

-- | Mahalanobis or ZCA whitening
zca :: U.Vector (Double, Double) -> U.Vector (Double, Double)
zca points = U.map (`dot'` w) $ points -- U.zip xs' ys'
  where
    w = t $ evec `dot` d `dot` t evec
    d = (1 / sqrt (eval1 + eps), 0, 0, 1 / sqrt (eval2 + eps))
    (evec, (eval1, eval2)) = eigendecomposition (sx, sxy, sxy, sy)
    eps = 10e-6
    sx = U.sum (U.zipWith (*) xs' xs') / n
    sy = U.sum (U.zipWith (*) ys' ys') / n
    sxy = U.sum (U.zipWith (*) xs' ys') / n
    xs' = U.map (\x -> x - mx) xs
    ys' = U.map (\x -> x - my) ys
    mx = mean xs
    my = mean ys
    n = fromIntegral $ U.length points - 1
    (xs, ys) = U.unzip points
{-# INLINE zca #-}

type Mat2x2 = (Double, Double, Double, Double)

-- | Eigendecomposition of 2 x 2 symetric matrix
eigendecomposition :: Mat2x2
                   -> ( Mat2x2   -- ^ Eigenvectors in the columns
                      , (Double, Double) ) -- ^ Eigenvalues
eigendecomposition (a,b,_,d) = 
    let (x1, y1) = evector eval1
        (x2, y2) = evector eval2
    in ((x1, x2, y1, y2), (eval1, eval2))
  where
    eval1 = (a + d + s) / 2
    eval2 = (a + d - s) / 2
    s = sqrt $ (a - d)**2 + 4 * b**2
    evector lambda = 
        let denominator = sqrt $ b**2 + (lambda - d)**2
        in ((lambda - d) / denominator, b / denominator)
{-# INLINE eigendecomposition #-}

dot :: Mat2x2 -> Mat2x2 -> Mat2x2
dot (a,b,c,d) (a',b',c',d') =
    ( a * a' + b * c'
    , a * b' + b * d'
    , c * a' + d * c'
    , c * b' + d * d' )
{-# INLINE dot #-}

dot' :: (Double, Double) -> Mat2x2 -> (Double, Double)
dot' (x,y) (a,b,c,d) = (x * a + y * c, x * b + y * d)
{-# INLINE dot' #-}

t :: Mat2x2 -> Mat2x2
t (a,b,c,d) = (a,c,b,d)
{-# INLINE t #-}
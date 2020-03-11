{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE TypeFamilies #-}
module Bio.ChromVAR
    ( computeDeviation
    , getBackgroundPeaks
    ) where

import qualified Eigen.Matrix as E
import qualified Eigen.SparseMatrix as ES
import qualified Eigen.Arithmetic as A
import Eigen.Matrix (Matrix)
import Eigen.SparseMatrix (SparseMatrix)
import Data.Singletons
import Data.Singletons.TypeLits
import Data.List
import Data.Maybe
import Conduit

import Bio.ChromVAR.Background

computeDeviation :: forall io. Monad io
                 => Int -- ^ The number of motifs
                 -> Int -- ^ The number of peaks
                 -> [Double]  -- ^ Expectation
                 -> [[Int]]   -- ^ Background assignment
                 -> [(Int, Int, Double)]   -- ^ peak by motif
                 -> ConduitT (Int, [(Int, Int, Double)]) ([[Double]], [[Double]]) io ()
computeDeviation nMotif nPeak expectation bg motifs =
    withSomeSing (fromIntegral nMotif) $ \(SNat :: Sing motif) -> 
        withSomeSing (fromIntegral nPeak) $ \(SNat :: Sing peak) -> 
            let e = fromJust $ E.fromList [expectation]
                background = flip map bg $ \b -> ES.transpose $ ES.fromList $ zipWith (\i j -> (i, j-1, 1)) [0..]  b
                peakByMotif = ES.fromList motifs :: SparseMatrix peak motif Double
                f (n, cell_by_peak) = withSomeSing (fromIntegral n) $ \(SNat :: Sing n) ->
                    let p = ES.fromList cell_by_peak :: SparseMatrix n peak Double
                        (mean, sd) = computeDeviation' e background peakByMotif p
                    in (E.toList mean, E.toList sd)
            in mapC f

-- | Compute the Bias-corrected deviations and z-scores.
computeDeviation' :: forall n m r . (SingI n, SingI m, SingI r)
                  => Matrix 1 m Double
                  -> [SparseMatrix m m Double]   -- ^ background assignment
                  -> SparseMatrix m r Double  -- ^ peak by motifs
                  -> SparseMatrix n m Double   -- cell by peak
                  -> (Matrix n r Double, Matrix n r Double)
computeDeviation' expectation bk motif_by_peak cell_by_peak = (deviations, deviations `div'` sd)
  where
    deviations = rawDev `E.sub` mean
    rawDev = (ES.toMatrix observed `E.sub` expected) `div'` expected
    observed = cell_by_peak `A.mul` motif_by_peak
    expected = readCount `A.mul` expectation `A.mul` motif_by_peak
    readCount = cell_by_peak `A.mul` E.ones
    mean = E.map (/ (fromIntegral $ length xs)) $ foldl1' E.add xs
    sd = E.map (sqrt . (/ (fromIntegral $ length xs - 1))) $ foldl1' E.add $ map (\x -> E.map (**2) $ x `E.sub` mean) xs
    xs = flip map bk $ \b -> 
        let sampled = cell_by_peak `A.mul` b `A.mul` motif_by_peak
            sampled_expected = readCount `A.mul` expectation `A.mul` b `A.mul` motif_by_peak
        in (ES.toMatrix sampled `E.sub` sampled_expected) `div'` sampled_expected
{-# INLINE computeDeviation' #-}

div' :: (SingI n, SingI m)
     => Matrix n m Double
     -> Matrix n m Double
     -> Matrix n m Double
div' m1 m2 = fromJust $ E.fromList $ (zipWith (zipWith (/))) (E.toList m1) $ E.toList m2
{-# INLINE div' #-}
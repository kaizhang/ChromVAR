{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE TypeFamilies #-}
module Bio.ChromVAR
    ( computeDeviation
    , computeDeviation'
    , getBackgroundPeaks
    ) where

import Data.Singletons
import Data.Singletons.TypeLits
import Data.List
import Conduit
import qualified Data.Vector.Storable as VS
import qualified Data.Matrix.Static.Dense as D
import qualified Data.Matrix.Static.Sparse as S
import Data.Matrix.Static.LinearAlgebra

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
            let e = D.matrix [expectation]
                background = flip map bg $ \b -> D.transpose $ S.fromTriplet $ zipWith (\i j -> (i, j, 1)) [0..]  b
                peakByMotif = S.fromTriplet motifs :: SparseMatrix peak motif Double
                f (n, cell_by_peak) = withSomeSing (fromIntegral n) $ \(SNat :: Sing n) ->
                    let p = S.fromTriplet cell_by_peak :: SparseMatrix n peak Double
                        (mean, sd) = computeDeviation' e background peakByMotif p
                    in (map VS.toList $ D.toRows mean, map VS.toList $ D.toRows sd)
            in mapC f

-- | Compute the Bias-corrected deviations and z-scores.
computeDeviation' :: forall n m r . (SingI n, SingI m, SingI r)
                  => Matrix 1 m Double
                  -> [SparseMatrix m m Double]   -- ^ background assignment
                  -> SparseMatrix m r Double  -- ^ peak by motifs
                  -> SparseMatrix n m Double   -- cell by peak
                  -> (Matrix n r Double, Matrix n r Double)
computeDeviation' expectation bk motif_by_peak cell_by_peak = (deviations, deviations / sd)
  where
    deviations = rawDev - mean
    rawDev = (D.convertAny observed - expected) / expected
    observed = cell_by_peak %*% motif_by_peak
    expected = readCount %*% expectation %*% motif_by_peak
    readCount = cell_by_peak %*% D.replicate 1
    mean = D.map (/ (fromIntegral $ length xs)) $ foldl1' (+) xs
    sd = D.map (sqrt . (/ (fromIntegral $ length xs - 1))) $
        foldl1' (+) $ map (\x -> D.map (**2) $ x - mean) xs
    xs = flip map bk $ \b -> 
        let sampled = cell_by_peak %*% b %*% motif_by_peak
            sampled_expected = readCount %*% expectation %*% b %*% motif_by_peak
        in (D.convertAny sampled - sampled_expected) / sampled_expected
{-# INLINE computeDeviation' #-}
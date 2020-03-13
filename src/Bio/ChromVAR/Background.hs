{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables   #-}
module Bio.ChromVAR.Background (getBackgroundPeaks) where

import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import qualified Data.Matrix.Unboxed as MU
import qualified Data.Vector.Unboxed.Mutable as UM
import Control.Monad
import Data.Ord
import Data.List
import Control.Arrow
import qualified Eigen.Matrix as E

import Statistics.Sample
import System.Random.MWC
import System.Random.MWC.Distributions (categorical)

import Bio.ChromVAR.Utils

type PeakGroup = ( MU.Matrix Double   -- ^ Group to group distance
                 , V.Vector (U.Vector Int)  -- ^ Peaks in each group
                 , U.Vector Int       -- ^ Peak-group membership
                 )

getBackgroundPeaks :: Int -> U.Vector (Double, Double) -> GenIO -> IO [U.Vector Int]
getBackgroundPeaks n xs gen =
    let pg = mkPeakGroup xs
    in replicateM n $ getBackgroundPeak gen pg

getBackgroundPeak :: GenIO -> PeakGroup -> IO (U.Vector Int)
getBackgroundPeak gen (weightMat, bins, membership) = U.generateM n $ \i -> do
    let ws = weightMat `MU.takeRow` (membership U.! i)
    grp <- categorical ws gen
    let peaks = bins V.! grp
    idx <- uniformR (0, U.length peaks - 1) gen
    return $ peaks U.! idx
  where
    n = U.length membership
{-# INLINE getBackgroundPeak #-}

-- | Make peak groups.
mkPeakGroup :: U.Vector (Double, Double) -> PeakGroup
mkPeakGroup raw = (weights, V.fromList groups, membership)
  where
    transformed = E.withList (map (\(x,y) -> [x,y]) $ U.toList raw) $ \mat ->
        U.fromList $ map (\[x,y] -> (x,y)) $ E.toList $ whiten Cholesky mat
    weights = MU.generate (U.length points, U.length points) $ \(i,j) ->
        weight (points U.! i) (points U.! j)
      where
        points = U.fromList $ flip map groups $ \is -> mean *** mean $
            U.unzip $ U.map (transformed U.!) is
    membership = U.create $ do
        v <- UM.new $ U.length transformed
        forM_ (zip [0..] groups) $ \(x, is) -> U.forM_ is $ \i ->
            UM.unsafeWrite v i x
        return v
    groups = go [] 0 0 $ sortBy (comparing fst) $ zip (U.toList transformed) [0..]
      where
        go acc i j (((x,y), idx) : rest)
            | null acc || i' == i || j' == j = go (idx : acc) i' j' rest
            | i' > i || j' > j = U.fromList acc : go [idx] i' j' rest
            | otherwise = error "Impossible"
          where
            i' = truncate $ (x - x_min) / x_step :: Int
            j' = truncate $ (y - y_min) / y_step :: Int
        go _ _ _ _ = []
    (xs, ys) = U.unzip transformed
    x_min = U.minimum xs
    x_max = U.maximum xs
    x_step = (x_max - x_min) / n
    y_min = U.minimum ys
    y_max = U.maximum ys
    y_step = (y_max - y_min) / n
    n = 50
{-# INLINE mkPeakGroup #-}
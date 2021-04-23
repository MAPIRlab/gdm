import numpy as np
from scipy.stats import multivariate_normal
from gdm.common.gdm import NormalGasDistributionMapper
from gdm.common import DiscreteScalarMap, Observation,  CellConnectivityMap,Lattice2DScalar

##============================================================================

DEFAULT_KERNEL_STD = 0.38
DEFAULT_SCALING_STD = 4.956
DEFAULT_WIND_STRECH = 0.33

##============================================================================

def windCovariance(std, w_speed, w_dir, strech):
	assert(std>0)
	a = std + strech * w_speed
	b = (std * std) / (std + strech * w_speed)
	E = np.array(((a * a, 0), (0, b * b)))

	if (w_speed > 0.01):
		R = np.array(((np.cos(w_dir), -np.sin(w_dir)), (np.sin(w_dir), np.cos(w_dir))))
		R_inv = np.linalg.inv(R)
		E_inv = np.linalg.inv(E)
		ER_inv = R_inv.dot(E_inv.dot(R))
		ER = np.linalg.inv(ER_inv)
		E = ER

	assert E[0,0] >= 0
	assert E[1, 1] >= 0

	return E



##============================================================================
class KDM_VW(NormalGasDistributionMapper):

	## CONSTRUCTORS ------------------------------------------------------------

	def __init__(self,
				 map,
	             kernel_std=DEFAULT_KERNEL_STD,
				 scaling_std=DEFAULT_SCALING_STD,
				 wind_strech=DEFAULT_WIND_STRECH):

		## Parameter check
		assert(kernel_std > 0)
		assert(scaling_std > 0)
		assert(wind_strech > 0)

		## Init base
		NormalGasDistributionMapper.__init__(self, size=map.size, resolution=map.resolution)

		## Member variables
		self._k_s = kernel_std
		self._o_s = scaling_std
		self._strech = wind_strech
		self._boundary = 10  # This sets a boundary around the point at which we are working to avoid updating the whole map. Smaller = faster computaiton

		self._omega = DiscreteScalarMap(size=self.size, resolution=self.resolution)
		self._R = DiscreteScalarMap(size=self.size, resolution=self.resolution)
		self._alpha = DiscreteScalarMap(size=self.size, resolution=self.resolution)

		self._variance = DiscreteScalarMap(size=self.size, resolution=self.resolution)
		self._sigma_gas = DiscreteScalarMap(size=self.size, resolution=self.resolution)
		self._variance = DiscreteScalarMap(size=self.size, resolution=self.resolution)
		self._variance_corrected = DiscreteScalarMap(size=self.size, resolution=self.resolution)
		self._max_omega = DiscreteScalarMap(size=self.size, resolution=self.resolution)
		self._variance_computed = 0

		return



	## PRIVATE -----------------------------------------------------------------




	## BASE CLASS IMPLEMENTATION -----------------------------------------------

	def _estimate(self):

		## INITI VARIABLES
		omega = np.zeros(self.shape)
		R = np.zeros(self.shape)
		total_gas = 0
		num_samples = 0

		## FOR EACH SAMPLE
		for sample in self._observations:
			total_gas += sample.gas
			num_samples += 1

			## COMPUTE COVARIANCE MATRIX DEPENDING ON WIND
			wind_speed = np.sqrt(sample.wind[0] ** 2 + sample.wind[1] ** 2)
			wind_direction = np.arcsin(sample.wind[0] / (wind_speed + 0.0001))
			covariance = windCovariance(self._k_s, wind_speed, wind_direction, self._strech)

			## COMPUTE REGION OF MAP WE WANT TO UPDATE (around sample)
			sample_cell = self._convertPositionToCell(sample.position)
			p = self._boundary
			i_scale = (p * max(np.sqrt(abs(covariance[0, 0])), np.sqrt(abs(covariance[1, 0])))) / self._k_s
			j_scale = (p * max(np.sqrt(abs(covariance[0, 1])), np.sqrt(abs(covariance[1, 1])))) / self._k_s
			start_i = max(0, int(sample_cell[0] - i_scale * self._k_s / self.resolution))
			start_j = max(0, int(sample_cell[1] - j_scale * self._k_s / self.resolution))
			end_i = min(self.shape[0], int(sample_cell[0] + 1 + i_scale * self._k_s / self.resolution) + 1)
			end_j = min(self.shape[1], int(sample_cell[1] + 1 + j_scale * self._k_s / self.resolution) + 1)

			num_cells_i = end_i - start_i
			num_cells_j = end_j - start_j
			index_i = [i for i in range(start_i, end_i)]
			index_j = [j for j in range(start_j, end_j)]
			cells = np.array(np.meshgrid(index_i, index_j)).T.reshape(-1, 2)
			distances = self.resolution * (cells - sample_cell)

			f = multivariate_normal.pdf(distances, cov=covariance).reshape(num_cells_i, num_cells_j)
			omega[start_i:end_i, start_j:end_j] += f
			R[start_i:end_i, start_j:end_j] += f * sample.gas


		## COMPUTE
		alpha = 1 - np.exp(-(omega ** 2) / (self._o_s ** 2))
		avg_gas = total_gas / (num_samples + 0.0001)
		gas = alpha * (R / (omega + 0.0001)) + (1 - alpha) * avg_gas

		## STORE DATA
		self._gas.loadMatrix(gas)
		self._R.loadMatrix(R)
		self._omega.toMatrix()
		self._omega.loadMatrix(omega)
		self._alpha.loadMatrix(alpha)		
		self._variance_computed = 0

		return self


	def _getCell(self, cell):
		gas = self.getGasEstimate().getCell(cell)
		position = self._convertCellToPosition(cell)
		o = Observation(position=position, gas=max(0.0, gas), data_type='gas')
		return o



	def _setCell(self, cell, value):
		assert False



	def plot(self):
		self._gas.plot()
		assert False


	def _normalize(self):
		assert False
		#return self


	def _triggerUpdateObstacleMap(self):
		return self


	def _computeUncertainty(self):
		a = self._alpha.toMatrix()
		u = 1/(a+0.00001) - 1
		u[u<0]=0.00001

		assert not np.isnan(u).any()
		assert u.min() > 0, u.min()
		assert u.max() < 99999999, u.max()

		self._gas_uncertainty.loadMatrix(u)
		return self


	def _getGasUncertainty(self):
		return self._gas_uncertainty



if __name__ == "__main__":
	from gdm.common.environments.corridor_1 import corridor_1

	kdmvw = KDM_VW(corridor_1)

	## DEFINE EXPLORATION PATH
	path = []
	path_checkpoints = ((0.5, 3.5), (14.5, 3.5), (14.5, 4.5), (0.5, 4.5))
	for i in range(0, len(path_checkpoints) - 1):
		checkpoint_1 = path_checkpoints[i]
		checkpoint_2 = path_checkpoints[i + 1]
		path += corridor_1.getStraightPathBetweenPositions(checkpoint_1, checkpoint_2)

	#path = [(2,4.5), (5,4.5), (2.5, 2.5)]

	## GET SAMPLES
	observations = corridor_1.getObservation(path)

	kdmvw.addObservation(observations)
	kdmvw.estimate()
	kdmvw.getGasEstimate().plot(vmax=1)
	kdmvw.getGasUncertainty().plot(vmax=2.5)


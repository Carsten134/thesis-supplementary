from torch.fft import fft2, fftshift
from torch import normal, tensor
import torch
import torch.nn as nn
import torch.nn.functional as F
import matplotlib.pyplot as plt


def randomize(x, y, times):
  # generate permutations
  perm = torch.bernoulli(torch.full((times, *x.shape), .5))
  perm_n = (perm == 0).to(torch.float16)

  # generate Tensor of shape (times, xH, xW)
  x_copy = x.repeat(times, 1, 1)
  y_copy = y.repeat(times, 1, 1)

  # do permutations all at once
  x_rand = x_copy * perm + y_copy * perm_n
  y_rand = x_copy * perm_n + y_copy * perm
  
  return x_rand.reshape(times, 1, *x.shape), y_rand.reshape(times, 1, *y.shape)  

def fourier_freq(N):
  return 2 * torch.pi * torch.arange(-(N - 1) // 2, N // 2) / N

def I(x):
  res = abs(fft2(x))**2 / (x.shape[0] * x.shape[1] * torch.var(x))
  return fftshift(res)

def k_2d_bp(N, M, hr = .2, hc = .2):
  C = 8/(torch.lgamma(tensor([.5])).exp()**2 * 4 * torch.pi**2 * hr * hc)
  omega_N = (fourier_freq(N) / hr)**2
  omega_M = (fourier_freq(M) / hc)**2
  
  omega_N = omega_N[omega_N < torch.pi**2]
  omega_M = omega_M[omega_M < torch.pi**2]

  omega_N = omega_N.reshape(omega_N.shape[0], 1)
  omega_M = omega_M.reshape(1, omega_M.shape[0])

  theta_2 = omega_N @ torch.ones(omega_M.shape[0]).reshape(1, omega_M.shape[0]) \
    + torch.ones(omega_N.shape[1]).reshape(omega_N.shape[1], 1) @ omega_M
  theta_2[theta_2 >= torch.pi**2] = torch.pi**2
  return C*(1-theta_2/torch.pi**2)


class BPKernel(nn.Module):
  def __init__(self, N, M, hr=.2, hc=.2, stride = 1):
      super().__init__()
      # build 2D kernel (numpy or torch) of shape (kH,kW)
      K = k_2d_bp(N, M, hr, hc)                  
      K = torch.as_tensor(K, dtype=torch.float32)

      # register it as a buffer of shape (1,1,kH,kW)
      self.register_buffer('K', K.view(1,1,*K.shape))
      # compute padding once
      self.padding = K.shape[0] // 2
      self.stride = stride

  def forward(self, x):
      # x: (batch,1,H,W); same device as self.K
      return F.conv2d(x, self.K, padding=self.padding, stride = self.stride) / (x.shape[2] * x.shape[3])


class gridMA(nn.Module):
  def __init__(self, K: tensor):
    super().__init__()
    # expecting square kernel of shape (KH, KH)
    K_buffer = torch.as_tensor(K, dtype= torch.float32)

    # register K in the buffer, such that it can be moved to a CUDA device
    self.register_buffer("K", K_buffer.reshape(1,1,*K_buffer.shape))
    self.padding = K.shape[0] // 2
    
  def forward(self, N, M, T = 1):
    N_tilde = N + 2 * self.padding
    M_tilde = M + 2 * self.padding
    eps = normal(0, 1, size = (T, N_tilde, M_tilde))
    eps = eps.reshape(T,1, N_tilde, M_tilde)
    return F.conv2d(eps, self.K)

def wold_w_q_fast(N, phi_1, phi_2, phi_3):
    wold_weights = torch.full((N, N),0, dtype=torch.float32)
    wold_weights[:, N -1] = phi_3 **torch.arange(N)
    wold_weights[0, [i for i in reversed(range(N))]] = phi_1**torch.arange(N)

    phi = torch.tensor([phi_1, phi_2, phi_3]).reshape(1, 3)
    for i in range(1, N):
        for j in reversed(range(N-1)):
            wold_weights[i, j] = phi @ torch.tensor([wold_weights[i,j+1], wold_weights[i-1, j+1], wold_weights[i-1, j]]).reshape(3, 1)
    
    # position wold weights in the left corner

    res = torch.zeros((2*N -1, 2*N - 1))
    res[(N-1):(2*N-1), 0:N] = wold_weights
    return res

def wold_w_h_fast(N, phi_1, phi_2, phi_3, phi_4):
  phi = torch.tensor([phi_1, phi_2, phi_3, phi_4]).reshape(1,4)

  # diagonal and horizontal weights
  res = torch.full((2*N - 1, 2*N - 1),0, dtype=torch.float32)
  res[[i for i in range(N-1, 2*N-1)], [i for i in range(N-1, 2*N - 1)]] = phi_4 **torch.arange(N)
  res[N-1, [i for i in reversed(range(N))]] = phi_1**torch.arange(N)

  for i in range(N, 2*N-1):
    for j in reversed(range(0, i)):
      if i == 0:
        res = [i, j] = phi @ torch.tensor([res[i,j+1], res[i-1, j+1], res[i-1, j], 0]).reshape(4, 1)
      else:
        res[i, j] = phi @ torch.tensor([res[i,j+1], res[i-1, j+1], res[i-1, j], res[i-1, j-1]]).reshape(4, 1)
  return res


# now implementing the quarter plane sample
class PlaneSampler(torch.nn.Module):
  def __init__(self, wold_size, phi_1, phi_2, phi_3, phi_4 = 0):
    super().__init__()
    # constructing wold weights 
    K_buffer = wold_w_q_fast(wold_size, phi_1,  phi_2, phi_3) if phi_4 == 0 else wold_w_h_fast(wold_size, phi_1, phi_2, phi_3, phi_4)
    
    # register K in the buffer, such that it can be moved to a CUDA device
    self.register_buffer("K", K_buffer.reshape(1,1,*K_buffer.shape))
    self.padding = K_buffer.shape[0] // 2

  def forward(self, N, M, T = 1):
    N_tilde = N + 2 * self.padding
    M_tilde = M + 2 * self.padding
    eps = normal(0, 1, size = (T, N_tilde, M_tilde))
    eps = eps.reshape(T,1, N_tilde, M_tilde)
    return F.conv2d(eps, self.K)
  
  def plot(self, N = 100):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    x = self.forward(N, N)[0][0]

    fig.dpi = 220
    ax1.imshow(x)
    ax1.set_title("Sample")
    ax2.imshow(I(x))
    ax2.set_title("Spectrum")
    ax3.imshow(self.K[0][0])
    ax3.set_title("Wold weights")

    plt.show()


# now implementing the test
def phi_n_star(x, y, B, alpha, hr = .2, hc = .2, stride = 1):

  # instanciate a kernel 
  bpkernel = BPKernel(*x.shape, hr, hc, stride).to(torch.device("cuda") if torch.cuda.is_available() else torch.float32)
  # first compute the periodogramms
  I_x = I(x).reshape(1,1, *x.shape).to(torch.device("cuda") if torch.cuda.is_available() else torch.float32)
  I_y = I(y).reshape(1,1, *x.shape).to(torch.device("cuda") if torch.cuda.is_available() else torch.float32)
  I_tilde = .5 * (I_x + I_y).to(torch.device("cuda") if torch.cuda.is_available() else torch.float32)

  # secondly compute Tn
  diff_x = bpkernel(I_x - I_tilde)**2
  diff_y = bpkernel(I_y - I_tilde)**2
  Tn = diff_x.sum() + diff_y.sum()
  Tn = Tn * (2*torch.pi)**2 / (diff_x.shape[2] * diff_x.shape[3])

  # then compute Tn_star
  I_tilde = I_tilde.reshape(*x.shape)
  I_x_rand, I_y_rand = randomize(I_x.reshape(*x.shape), I_y.reshape(*x.shape), B)
  
  diff_rand_x = bpkernel(I_x_rand - I_tilde.repeat(B,1, 1, 1))**2
  diff_rand_y = bpkernel(I_y_rand - I_tilde.repeat(B,1, 1, 1))**2
  Tn_star = diff_rand_x.sum((2,3)) + diff_rand_y.sum((2, 3))
  Tn_star = Tn_star  * (2*torch.pi)**2 / (diff_rand_x.shape[2] * diff_rand_x.shape[3])
  
  return ((Tn_star < Tn).sum() / B > 1 - alpha).to(torch.float32), Tn, Tn_star
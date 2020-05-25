# Numerical_Toolkit_for_Approximating_Groundwater_Flow_Field_from_Scattered_Data

![Preview](https://numericalenvironmental.files.wordpress.com/2020/05/script-particle-tracks.jpg?w=816)

This Julia script imports and interpolates scattered hydraulic head and hydraulic conductivity data (for a single-layer aquifer) using ordinary kriging. The interpolation function serves as a basis for several groundwater flow assessment tools, including (1) a gridding function for subsequent plotting, (2) groundwater velocity vector estimation as a function of location, (3) forward and backward particle tracking, and (4) numerical integration of groundwater fluxes along a user-supplied polyline or polygon. Additional information and a summary of an example application can be found on my blog, https://numericalenvironmental.wordpress.com/2020/05/25/a-numerical-toolkit-for-summarizing-groundwater-flow-patterns-in-a-2-d-aquifer-from-scattered-data/.
The script runs under Julia 1.0 or later, and requires the following packages: DelimitedFiles, Statistics, Polynomials, and QuadGK. Text-based input files include the following:
* grid.txt – model boundaries (and grid resolution to support plotting)
* heads.txt – groundwater elevations over (x, y)
* logKs.txt – log hydraulic conductivity values over (x, y)
* particleStarts.txt – starting locations and migration direction for particle tracks; read if trackFlag is set = true
* trackParams.txt – porosity and particle step size; read if trackFlag is set = true
* polygon.txt – apex points for polyline/polygon for flux integration; read if polyFlag is set = true
* wells.txt – well locations and pumping rates; if wellFlag is set = true (pumping rates not used in current version)

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


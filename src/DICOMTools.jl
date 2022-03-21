
module DICOMTools

	export get_tag
	export get_transformation_matrix
	export dcm2array
	export load_dcms, load_dcms!
	export merge_multiseries
	export DICOMDict
	export DICOMSeries

	using DICOM

	# Sorry for the bad commenting/documenting ... when I find the time I'll add it

	# Define a types to store DICOMs found in a directory structure
	struct DICOMSeries
		paths::Vector{String} # Individual dicom files
		dcms::Vector{DICOM.DICOMData} # The read in dicom structures
		pos::Vector{Int64} # Position in the files where reading stopped (see enum WhatToRead)
	end
	const DICOMDict = Dict{String, DICOMSeries} # SeriesInstanceUID => DICOMSeries

	# This enum can be used to define up to which group to read
	@enum WhatToRead read_all read_uid read_voxelsize skip_pixels
	# Convert above enum to a DICOM group
	function get_max_group(what::WhatToRead)::UInt16
		# Get DICOM.dcm_parse's max_group argument to stop reading early
		max_group = (what == read_all) ? 0xffff :
					(what == read_uid) ? 0x0020 :
					(what == read_voxelsize) ? 0x0028 : 0x7FDF # last option must be skip_pixels
		return max_group
	end



	# Load Dicoms recursively
	function load_dcms(directory::String, what::WhatToRead=read_all)::DICOMDict

		max_group = get_max_group(what)

		dcmdict = DICOMDict()
		for (root, _, files) in walkdir(directory)
			for filename in files
				fullpath = joinpath(root, filename)
				DICOM.isdicom(fullpath) || continue
				local dcm, pos
				open(fullpath, "r") do file
					dcm = dcm_parse(file; max_group)
					if eof(file)
						pos = -1
					else
						pos = position(file) - sizeof(UInt16) # Group tag was read in, revert position
					end
				end
				uid = dcm[tag"SeriesInstanceUID"]
				if haskey(dcmdict, uid)
					push!(dcmdict[uid].paths, fullpath)
					push!(dcmdict[uid].dcms, dcm)
					push!(dcmdict[uid].pos, pos)
				else
					dcmdict[uid] = DICOMSeries([fullpath], [dcm], [pos])
				end
			end
		end
		return dcmdict
	end

	# Read tags retrospectively
	function load_dcms!(dcmdict::DICOMDict, what::WhatToRead=read_all)::DICOMDict
		max_group = get_max_group(what)
		for series in values(dcmdict)
			load_dcms!(series, max_group)
		end
		return dcmdict
	end
	function load_dcms!(series::DICOMSeries, what::WhatToRead=read_all)::DICOMSeries
		max_group = get_max_group(what)
		load_dcms!(series::DICOMSeries, max_group)
		return series
	end
	function load_dcms!(series::DICOMSeries, max_group::UInt16)::DICOMSeries
		for (path, dcm, pos) in zip(series.paths, series.dcms, series.pos)
			pos == -1 && continue
			open(path, "r") do file
				seek(file, pos)
				DICOM.read_body!(file, dcm; max_group)
			end
		end
		return series
	end


	function filter(f::Function, dcmdict::DICOMDict)
		Base.filter(series -> f(series.second.dcms), dcmdict)
	end
	function filter!(f::Function, dcmdict::DICOMDict)
		Base.filter!(series -> f(series.second.dcms), dcmdict)
	end


	function get_tag(
		multiseries::AbstractVector{<: DICOMSeries},
		tag::Tuple{UInt16,UInt16},
		outtype::Val{T} = Val(typeof(multiseries[1].dcms[1][tag]))
	) where T
		T[ series.dcms[1][tag] for series in multiseries ]
	end
	function get_tags(
		multiseries::AbstractVector{<: DICOMSeries},
		tags::NTuple{N, Tuple{UInt16,UInt16}},
		outtypes::Val{T} = Val(ntuple(i -> typeof(multiseries[1].dcms[1][tags[i]]), N))
	) where {N, T}
		T[ series.dcms[1][tag] for tag in tags for series in multiseries ]
	end


	function get_transformation_matrix(series::DICOMSeries)::Matrix{Float64}
		# This breaks if InstanceNumber tag does not equal slice number
		# Gaps between slices are handled correclty if volume is rectangular

		dcms = series.dcms

		# Get origin of first and second slice
		# Translation of first slice is origin of the volume
		# Translation of second minus translation of first is vector pointing in slice direction.
		local translation::Vector{Float64}, slice_vector::Vector{Float64} 
		let first_slice = 0, second_slice = 0
			for (i, dcm) in enumerate(dcms)
				instance_number = dcm[(0x0020, 0x0013)]
				if instance_number == 1
					first_slice = i
				elseif instance_number == 2
					second_slice = i
				end
				first_slice != 0 && second_slice != 0 && break
			end
			translation = dcms[first_slice][(0x0020, 0x0032)]
			slice_vector = dcms[second_slice][(0x0020, 0x0032)] - translation
			# Do not normalise, so that slice thickness is included
		end

		# Get "step"-vectors in row and column direction
		local row_vector::Vector{Float64}, column_vector::Vector{Float64}
		let orientation = dcms[1][(0x0020, 0x0037)], pixelsize = dcms[1][(0x0028, 0x0030)]
			row_vector = pixelsize[1] .* orientation[1:3]
			column_vector = pixelsize[2] .* orientation[4:6]
		end

		transformation = Matrix{Float64}(undef, 4, 4)
		transformation[1:3, 1] = row_vector
		transformation[1:3, 2] = column_vector
		transformation[1:3, 3] = slice_vector
		transformation[1:3, 4] = translation
		transformation[4, 1:3] .= 0
		transformation[4, 4] = 1
		return transformation
	end



	function dcm2array(
		series::DICOMSeries,
		outtype::Val{T} = Val(Float64)
	)::Array{T, 3} where T
		# Note: This breaks if InstanceNumber tag does not equal slice number

		dcms = series.dcms

		# Get size
		reference = dcms[1]
		array_size = (
			get(reference, tag"Rows", 1),
			get(reference, tag"Columns", 1),
			maximum( get(dcm, tag"InstanceNumber", 1) for dcm in dcms )
		)

		# Fill array
		voxels = Array{T, 3}(undef, array_size...)
		for dcm in dcms
			data = dcm[tag"PixelData"]
			voxels[:, :, dcm[tag"InstanceNumber"]] = data
		end

		return voxels
	end






	"""

		Get an array extracted from multiple DICOM series.
		For example useful for MRI: fitting relaxation curves to series
		acquired with different inversion/echo times.

	"""
	function merge_multiseries(
		multiseries::AbstractVector{<: DICOMSeries},
		data_type::Val{T} = Val(eltype(multiseries[1].dcms[1]))
	)::Array where T
		# TODO: Would expect to return Array{..., 4} because dcm2array gives 3D volume.
		# However what happens with timeseries in 3D volumes?

		# Get arrays from dcms
		tojoin = [ dcm2array(dcms) for dcms in multiseries ]

		# Join arrays
		# Rearrange all volumes into a single array, first axis iterates over parameters,
		# this makes access to signal(parameter) easier.
		# Create empty array from the size of the first element in tojoin
		arrays = Array{
			T,
			1 + ndims(tojoin[1]) # parameter (1) and spatial dimensions
		}(
			undef,
			length(tojoin),
			size(tojoin[1])... # Use first element as reference
		)
		array_size = tuple((size(arrays, d) for d = 2:4)...)
		for (i, array) in enumerate(tojoin)
			@assert size(array) == array_size # Does size match?
			arrays[i, :, :, :] = array
		end
		return arrays
	end

end


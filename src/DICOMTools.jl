
module DICOMTools

	export gettag
	export dcm2array
	export load_dcms
	export load_pixeldata!
	export merge_multiseries

	using DICOM



	function gettag(
		dcmdict::AbstractVector{<: AbstractVector{DICOM.DICOMData}},
		tag::Tuple{UInt16,UInt16};
		outtype::Type = Nothing
	)
		# Get type
		if outtype <: Nothing
			outtype = eltype(dcmdict[1][1][tag])
		end

		# Get values
		return outtype[ dcms[1][tag] for dcms in dcmdict ]
	end



	function dcm2array(
		dcms::AbstractVector{<: DICOM.DICOMData};
		outtype::Type = Nothing
	)::Array{<: Real, 3}
		
		# Get size
		reference = dcms[1]
		array_size = collect(
			get(reference, tag, 1)
			for tag in (tag"Rows", tag"Columns", tag"NumberOfSlices")
		)

		# Get type
		if outtype <: Nothing
			outtype = eltype(reference[tag"PixelData"])
		end

		# Fill array
		voxels = Array{outtype, 3}(undef, array_size...)
		for dcm in dcms
			data = dcm[tag"PixelData"]
			voxels[:, :, dcm[tag"InstanceNumber"]] = data
		end

		return voxels
	end



	const DICOMDict = Dict{
		String,
		NamedTuple{
			(:paths, :dcms),
			Tuple{Vector{String}, Vector{DICOM.DICOMData}}
		}
	}

	function load_dcms(directory::String; skip_pixels::Bool=false)::DICOMDict
		# Should pixels be skipped?
		aux_vr = skip_pixels ? Dict((0x7FE0, 0x0010) => DICOM.empty_vr) : Dict()
		# TODO: Could use max_group, but then the question is what other tags the user might need.
		# Potentially it is good to just add aux_vr and max_group and hand them through to dcm_parse,
		# shifting the problem on the users' side

		dcmdict = DICOMDict()
		for (root, _, files) in walkdir(directory)
			for filename in files
				fullpath = joinpath(root, filename)
				dcm = dcm_parse(fullpath; aux_vr)
				uid = dcm[tag"SeriesInstanceUID"]
				if haskey(dcmdict, uid)
					push!(dcmdict[uid].paths, fullpath)
					push!(dcmdict[uid].dcms, dcm)
				else
					dcmdict[uid] = (paths=[fullpath], dcms=[dcm])
				end
			end
		end
		return dcmdict
	end



	function read_tag!(dcm::DICOM.DICOMData, tag::Tuple{UInt16, UInt16}, filename::String)
		# Remove tag from vr
		delete!(dcm.vr, tag)

		# Get starting position in file
		start = dcm[tag]

		# Open file and read
		open(filename, "r") do f
			seek(f, start)
			(gelt, data, vr) = DICOM.read_element(f, dcm)
			if gelt == DICOM.empty_tag
				error("Could not read element $tag")
			else
				dcm[gelt] = data
			end
		end
	end

	function load_pixeldata!(dcmdict::DICOMDict)
		for element in values(dcmdict)
			for (path, dcm) in zip(element.paths, element.dcms)
				read_tag!(dcm, tag"PixelData", path)
			end
		end
	end



	"""
		merge_multiseries(dcmdict, uids, tag, data_type=Nothing, tag_type=Nothing)

		Get an array extracted from multiple DICOM series.
		For example useful for MRI: fitting relaxation curves to series
		acquired with different inversion/echo times.

	"""
	function merge_multiseries(
		dcmdict::DICOMDict,
		uids::AbstractVector{<: String},
		tag::Tuple{UInt16, UInt16};
		data_type::Type = Nothing,
		tag_type::Type = Nothing
	)::Tuple{Array, Vector}
		# TODO: Add option to Register images?

		# Get the respective DICOMs
		selected = Vector{Vector{DICOM.DICOMData}}(undef, length(uids))
		for (i, uid) in enumerate(uids)
			element = get(dcmdict, uid, nothing)
			isnothing(element) && error("UID $uid not found in DICOM batch")
			selected[i] = element.dcms
		end

		# Extract tags
		tags = gettag(selected, tag; outtype=tag_type)

		# Get arrays from dcms
		tojoin = [ dcm2array(dcms) for dcms in selected ]

		# Join arrays
		# Rearrange all volumes into a single array, first axis iterates over parameters,
		# this makes access to signal(parameter) easier.
		# Create empty array from the size of the first element in tojoin
		if data_type <: Nothing
			data_type = eltype(tojoin[1])
		end
		arrays = Array{
			data_type,
			1 + ndims(tojoin[1]) # parameter (1) and spatial dimensions
		}(
			undef,
			length(tags),
			size(tojoin[1])... # Use first element as reference
		)
		array_size = tuple((size(arrays, d) for d = 2:4)...)
		for (i, array) in enumerate(tojoin)
			@assert size(array) == array_size # Does size match?
			arrays[i, :, :, :] = array
		end
		return arrays, tags
	end

end

